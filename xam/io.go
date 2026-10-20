package xam

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ── Constants ────────────────────────────────────────────────────────────────

const ChanCap = 5000
const IoBuf = 4 << 20 // 4 MiB

// ── Interfaces ────────────────────────────────────────────────────────────────

type AlignmentReader interface {
	Read() (*sam.Record, error)
	Header() *sam.Header
	Close() error
}

type AlignmentWriter interface {
	Write(*sam.Record) error
	Close() error
}

// ── Misc ────────────────────────────────────────────────────────────────

// Given a thread total, allots read and write threads, keeping read constant at 1
// and scaling write. Returns: read, write.
func IoThreads(threads int) (int, int) {
	return 1, max(1, threads-1)
}

// ── SAM wrappers (sam.Reader/Writer lack Close()) ─────────────────────────────

type SamReader struct {
	*sam.Reader
	closer io.Closer
}

func (s *SamReader) Close() error {
	return s.closer.Close()
}

type SamWriter struct {
	*sam.Writer
}

func (s *SamWriter) Close() error {
	return nil
}

// ── Helpers ───────────────────────────────────────────────────────────────────

func checkError(err error) {
	if err != nil {
		log.Fatalf("%v", err)
	}
}

func BamNotStdout(isSam bool) error {
	stat, err := os.Stdout.Stat()
	if err != nil {
		return err
	}
	if (stat.Mode()&os.ModeCharDevice) != 0 && !isSam {
		return fmt.Errorf("refusing to write BAM to the terminal, as this was likely unintended. Use -S or pipe to another tool/file")
	}
	return nil
}

// FileOrStdin resolves a positional argument to a file path or "-" for stdin.
func FileOrStdin(infile string) (string, error) {
	if len(infile) == 0 {
		stat, err := os.Stdin.Stat()
		if err != nil {
			return "", err
		}
		if (stat.Mode() & os.ModeCharDevice) != 0 {
			return "", fmt.Errorf("no input file given and stdin is not piped data")
		}
		return "-", nil
	}
	return infile, nil
}

// Create a formatted Djinn PG line to add to SAM header
func NewPG(hdr *sam.Header, cl string) *sam.Program {
	progs := hdr.Progs()

	referenced := make(map[string]bool, len(progs))
	seen := make(map[string]bool, len(progs))
	for _, p := range progs {
		if pp := p.Previous(); pp != "" {
			referenced[pp] = true
		}
		seen[p.UID()] = true
	}

	var prev string
	for _, p := range progs {
		if !referenced[p.UID()] {
			prev = p.UID()
			break
		}
	}

	uid := "djinn"
	for i := 1; seen[uid]; i++ {
		uid = fmt.Sprintf("djinn.%d", i)
	}

	return sam.NewProgram(uid, "djinn", cl, prev, "3.0")
}

// ── Reader channel ────────────────────────────────────────────────────────────

// openAlignmentReader opens inFile (or stdin if "-") and sniffs its magic
// bytes to decide between a BAM and SAM reader. isBAM tells the caller
// whether it, not SamReader.Close(), owns closing fh.
func openAlignmentReaderBuf(inFile string, buff, threads int) (r AlignmentReader, fh *os.File, bufReader *bufio.Reader, isBAM bool) {
	fh, err := os.Stdin, error(nil)
	if inFile != "-" {
		fh, err = os.Open(inFile)
		checkError(err)
	}

	bufReader = bufio.NewReaderSize(fh, buff)

	magic, err := bufReader.Peek(2)
	if err != nil {
		if fh != os.Stdin {
			fh.Close()
		}
		checkError(err)
	}

	isBAM = magic[0] == 0x1f && magic[1] == 0x8b
	if isBAM {
		br, err := bam.NewReader(bufReader, threads)
		checkError(err)
		r = br
	} else {
		sr, err := sam.NewReader(bufReader)
		checkError(err)
		r = &SamReader{Reader: sr, closer: fh}
	}
	return r, fh, bufReader, isBAM
}

// NewXamReaderChan opens a SAM or BAM file (or stdin if inFile == "-") and
// streams records into the returned channel. The returned AlignmentReader
// exposes the file header. The goroutine owns all cleanup.
func NewXamReaderChan(inFile string, cp, buff, threads int) (chan *sam.Record, AlignmentReader) {
	outChan := make(chan *sam.Record, cp)
	r, fh, _, isBAM := openAlignmentReaderBuf(inFile, buff, threads)

	go func() {
		// Only close fh here for BAM; SamReader.Close() owns it for SAM.
		if isBAM && fh != os.Stdin {
			defer fh.Close()
		}
		for {
			rec, err := r.Read()
			if err == io.EOF {
				close(outChan)
				return
			}
			if err != nil {
				close(outChan)
				checkError(err) // fatal — exits, no fallthrough
				return          // defensive
			}
			outChan <- rec
		}
	}()

	return outChan, r
}

// NewXamReaderChanBatch is the batched counterpart of NewXamReaderChan: the
// reader goroutine accumulates up to batchSize records before sending, so a
// consumer that just wants records in bulk (e.g. to re-chunk them) pays one
// channel send/receive per batch instead of one per record. Behaves
// identically otherwise, including a final short batch at EOF.
func NewXamReaderChanBatch(inFile string, cp, buff, threads, batchSize int) (chan []*sam.Record, AlignmentReader) {
	outChan := make(chan []*sam.Record, cp)
	r, fh, _, isBAM := openAlignmentReaderBuf(inFile, buff, threads)

	go func() {
		if isBAM && fh != os.Stdin {
			defer fh.Close()
		}
		batch := make([]*sam.Record, 0, batchSize)
		for {
			rec, err := r.Read()
			if err == io.EOF {
				if len(batch) > 0 {
					outChan <- batch
				}
				close(outChan)
				return
			}
			if err != nil {
				close(outChan)
				checkError(err)
				return
			}
			batch = append(batch, rec)
			if len(batch) >= batchSize {
				outChan <- batch
				batch = make([]*sam.Record, 0, batchSize)
			}
		}
	}()

	return outChan, r
}

// ── Writer channel ────────────────────────────────────────────────────────────

// openAlignmentWriter opens outFile (or stdout if "-") and builds a SAM or
// BAM writer over a buffered wrapper around it.
func openAlignmentWriter(outFile string, head *sam.Header, buff, threads int, uncompressed bool) (w AlignmentWriter, fh *os.File, bio *bufio.Writer) {
	fh, err := os.Stdout, error(nil)
	if outFile != "-" {
		fh, err = os.Create(outFile)
		checkError(err)
	}

	bio = bufio.NewWriterSize(fh, buff)

	if uncompressed {
		sw, err := sam.NewWriter(bio, head, sam.FlagDecimal)
		checkError(err)
		w = &SamWriter{Writer: sw}
	} else {
		bw, err := bam.NewWriterLevel(bio, head, 4, threads)
		checkError(err)
		w = bw
	}
	return w, fh, bio
}

// NewXamWriterChan writes records received on the returned channel to outFile
// (or stdout if outFile == "-"). Sends true on doneChan when the input channel
// is closed and all records have been flushed.
func NewXamWriterChan(outFile string, head *sam.Header, cp, buff, threads int, uncompressed bool) (chan *sam.Record, chan bool) {
	outChan := make(chan *sam.Record, cp)
	doneChan := make(chan bool)
	w, fh, bio := openAlignmentWriter(outFile, head, buff, threads, uncompressed)

	go func() {
		for rec := range outChan {
			checkError(w.Write(rec))
		}
		w.Close()
		checkError(bio.Flush()) // always flush
		if fh != os.Stdout {
			fh.Close()
		}
		doneChan <- true
	}()

	return outChan, doneChan
}

// NewXamWriterChanBatch is the batched counterpart of NewXamWriterChan: the
// caller sends slices of records instead of one record per channel op, so a
// producer that already has records in bulk (e.g. off a merge heap) pays one
// channel send per batch instead of one per record.
func NewXamWriterChanBatch(outFile string, head *sam.Header, cp, buff, threads int, uncompressed bool) (chan []*sam.Record, chan bool) {
	outChan := make(chan []*sam.Record, cp)
	doneChan := make(chan bool)
	w, fh, bio := openAlignmentWriter(outFile, head, buff, threads, uncompressed)

	go func() {
		for batch := range outChan {
			for _, rec := range batch {
				checkError(w.Write(rec))
			}
		}
		w.Close()
		checkError(bio.Flush())
		if fh != os.Stdout {
			fh.Close()
		}
		doneChan <- true
	}()

	return outChan, doneChan
}
