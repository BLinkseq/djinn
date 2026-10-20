// Package fqsort sorts paired-end FASTQ reads by linked-read barcode,
// mirroring the sort package's BAM/SAM k-way merge sort: read + pair,
// sort-and-spill chunks in parallel, then k-way merge the spilled chunks
// back into barcode order.
package sort

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"github.com/blinkseq/djinn/barcodes"
	"github.com/blinkseq/djinn/fastq"
	"io"
	"os"
	"slices"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

// mate is one read's sequence/quality plus its CASAVA marker verbatim
// (nil if the source record had none). The marker's own text matters,
// not just whether one was present: ToStandard and ToHaplotagging only
// ever read its first byte (the '1'/'2' mate digit), but ToStlfr,
// ToTellseq, and ToTenX all write the *whole* marker back out when it's
// longer than that single byte -- a new-style Illumina CASAVA field
// ("1:N:0:ATCGATCG") carries its own filter flag and index sequence that
// a synthesized placeholder can't reconstruct. Storing only a hasCasava
// bool here (an earlier version of this file did) silently replaced
// that real marker with a generic stand-in on round-trip for those
// three formats -- wrong for anything but ToStandard/ToHaplotagging.
type mate struct {
	seq, qual []byte
	casava    []byte
}

// fqPair is a self-contained (deep-copied) record: one mate for
// single-end input, two for paired -- every record within a single
// SortFqByBX call has the same number, fixed by how many files were
// given -- plus the barcode/validity resolved from whichever mate (or
// the lone mate, single-end) carried it. Unlike the xam sort's bxKeyOf,
// which can alias the BAM reader's own buffers because a record is
// re-encoded to its temp file immediately, fqPair must own its memory: a
// chunk of these sits in RAM across a sort before being spilled, long
// after fastx.Reader has reused its internal buffers for later records.
//
// id and comments are stored once, not once per mate: every format
// CoreFqParser detects carries its barcode/CASAVA marker as the only
// thing that structurally differs between mates, so once those are
// stripped out (which CoreFq already does), the remaining description
// text is identical across mates -- see resolveComments.
type fqPair struct {
	barcode  []byte // nil if no mate carried one
	id       []byte // shared base read ID (CASAVA suffix stripped)
	comments []byte // shared description text once CASAVA/barcode are stripped
	valid    bool
	mates    []mate // len 1 (single-end) or 2 (paired)
}

// pairKey returns p's sort key: its barcode, or missingBXKey when p has
// none (unset barcodes sort last, matching the xam sort's convention).
func pairKey(p *fqPair) []byte {
	if len(p.barcode) == 0 {
		return missingBXKey
	}
	return p.barcode
}

// resolveBarcode prefers mate 1's barcode, falling back to mate 2's; nil
// if neither mate's core parsed to one.
func resolveBarcode(c1 fastq.CoreFq, ok1 bool, c2 fastq.CoreFq, ok2 bool) []byte {
	if ok1 && len(c1.Barcode) > 0 {
		return c1.Barcode
	}
	if ok2 && len(c2.Barcode) > 0 {
		return c2.Barcode
	}
	return nil
}

// resolveValid mirrors resolveBarcode's preference, since a record with
// no parsed core has a meaningless Valid value.
func resolveValid(c1 fastq.CoreFq, ok1 bool, c2 fastq.CoreFq, ok2 bool) bool {
	if ok1 {
		return c1.Valid
	}
	if ok2 {
		return c2.Valid
	}
	return false
}

// resolveComments assumes both mates carry identical description text
// once their barcode/CASAVA markers are stripped -- true for every format
// CoreFqParser recognizes, since the barcode tag/marker is the only thing
// that differs structurally between mates. Prefers whichever side is
// non-empty, defaulting to mate 1, so that assumption's one failure mode
// (the two genuinely differ) still keeps *something* rather than
// silently losing both.
func resolveComments(c1, c2 fastq.CoreFq) []byte {
	if len(c1.Comments) > 0 {
		return c1.Comments
	}
	return c2.Comments
}

// casavaOne and casavaTwo are the literal CASAVA bytes ToStandard ever
// writes -- it only ever reads rec.CASAVA[0], the mate digit -- so a
// stored pair only needs to know whether each mate had a marker at all,
// not the marker's own text.
var casavaOne = []byte{'1'}
var casavaTwo = []byte{'2'}

// SortFqByBX reads fqs (one unpaired FASTQ file, or two sharing matching
// read IDs) as a record stream, sorts by linked-read barcode, and writes
// outPrefix+".R1.fq.gz" (and outPrefix+".R2.fq.gz" when fqs has two
// entries). The barcode format (haplotagging/stLFR/TELL-seq) is
// auto-detected from fqs[0], same as fastq.CoreFqParser's other callers.
// When none of those is detected and bcMap names a 10X barcode whitelist
// file, falls back to the 10X path instead -- same dispatch
// convert.ConvertFq uses: format detection first, an explicit whitelist
// as the fallback, since 10X's barcode lives inline in R1's sequence and
// has no marker a sampling pass could recognize on its own.
//
// Output is written in the same format it was read in (haplotagging stays
// haplotagging, etc.) unless standardize is true, in which case output
// uses fastq.CoreFq.ToStandard's generic BX:Z:/VX:i: form instead -- the
// same opt-in every other djinn fastq tool (count, extract, invalid,
// sample, unlinked) already leaves to the caller.
func SortFqByBX(fqs []string, outPrefix, tmpDir string, threads int, standardize bool, bcMap string) error {
	seq.ValidateSeq = false

	numMates := len(fqs)
	if numMates < 1 || numMates > 2 {
		return fmt.Errorf("SortFqByBX accepts one or two FASTQ files, got %d", numMates)
	}

	coreParser, fmtType, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return err
	}

	// parseMate1/parseMate2 are asymmetric only for 10X, where the
	// barcode lives solely in R1's sequence prefix: TenXR22Corefq never
	// finds one on its own, so its ok is always false, which is exactly
	// what resolveBarcode/resolveValid need to fall through to mate 1's
	// result below -- the same preference they already apply for every
	// other format's (rarer, but possible) single-mate-only barcode.
	parseMate1, parseMate2 := coreParser, coreParser
	if coreParser == nil {
		if bcMap == "" {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s. If this is 10X data, a barcode whitelist file must be provided to identify inline barcodes.", fqs[0])
		}
		bclist, err := barcodes.LoadTenXFromFile(bcMap)
		if err != nil {
			return err
		}
		parseMate1 = func(rec *fastx.Record) (fastq.CoreFq, bool) { return fastq.TenXR12Corefq(rec, bclist) }
		parseMate2 = func(rec *fastx.Record) (fastq.CoreFq, bool) { return fastq.TenXR22Corefq(rec), false }
		fmtType = "10x"
	}

	converter, err := formatConverter(fmtType, standardize)
	if err != nil {
		return err
	}

	r1Reader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
	}
	defer r1Reader.Close()

	var r2Reader *fastx.Reader
	if numMates == 2 {
		r2Reader, err = fastx.NewReader(seq.DNA, fqs[1], "")
		if err != nil {
			return fmt.Errorf("opening %s: %w", fqs[1], err)
		}
		defer r2Reader.Close()
	}

	if tmpDir != "" {
		if err := os.MkdirAll(tmpDir, 0o755); err != nil {
			return err
		}
	} else {
		tmpDir = "."
	}

	chunks := make(chan []fqPair, threads)
	results := make(chan string, threads)
	errCh := make(chan error, threads)

	var wg sync.WaitGroup
	for range threads {
		wg.Go(func() {
			for c := range chunks {
				path, err := fqsortAndSpill(c, tmpDir)
				if err != nil {
					select {
					case errCh <- err:
					default:
					}
					continue
				}
				results <- path
			}
		})
	}

	var tmpFiles []string
	var collectWG sync.WaitGroup
	collectWG.Go(func() {
		for path := range results {
			tmpFiles = append(tmpFiles, path)
		}
	})

	// producer: read record(s) on the calling goroutine (reading is
	// inherently serial here -- when paired, R1 and R2 must be advanced
	// in lockstep -- so there's no analogue to the xam sort's separate
	// reader channel).
	var readErr error
	var chunk []fqPair
	enqueue := func(p fqPair) {
		chunk = append(chunk, p)
		if len(chunk) >= chunkSize {
			chunks <- chunk
			chunk = nil
		}
	}

	if numMates == 1 {
		for {
			rec1, err1 := r1Reader.Read()
			if err1 == io.EOF {
				break
			}
			if err1 != nil {
				readErr = err1
				break
			}

			core1, ok1 := parseMate1(rec1)
			var bc []byte
			if ok1 {
				bc = core1.Barcode
			}
			enqueue(fqPair{
				barcode:  append([]byte(nil), bc...),
				id:       append([]byte(nil), core1.ID...),
				comments: append([]byte(nil), core1.Comments...),
				valid:    ok1 && core1.Valid,
				mates: []mate{{
					seq:    append([]byte(nil), core1.Seq...),
					qual:   append([]byte(nil), core1.Qual...),
					casava: append([]byte(nil), core1.CASAVA...),
				}},
			})
		}
	} else {
	readLoop:
		for {
			rec1, err1 := r1Reader.Read()
			rec2, err2 := r2Reader.Read()
			eof1, eof2 := err1 == io.EOF, err2 == io.EOF
			if eof1 && eof2 {
				break readLoop
			}
			if err1 != nil && !eof1 {
				readErr = err1
				break readLoop
			}
			if err2 != nil && !eof2 {
				readErr = err2
				break readLoop
			}
			if eof1 != eof2 {
				readErr = fmt.Errorf("%s and %s have different numbers of records (mate files out of sync)", fqs[0], fqs[1])
				break readLoop
			}

			core1, ok1 := parseMate1(rec1)
			core2, ok2 := parseMate2(rec2)
			if !bytes.Equal(core1.ID, core2.ID) {
				readErr = fmt.Errorf("mismatched read IDs between mates: %q (R1) vs %q (R2)", core1.ID, core2.ID)
				break readLoop
			}

			enqueue(fqPair{
				barcode:  append([]byte(nil), resolveBarcode(core1, ok1, core2, ok2)...),
				id:       append([]byte(nil), core1.ID...),
				comments: append([]byte(nil), resolveComments(core1, core2)...),
				valid:    resolveValid(core1, ok1, core2, ok2),
				mates: []mate{
					{
						seq:    append([]byte(nil), core1.Seq...),
						qual:   append([]byte(nil), core1.Qual...),
						casava: append([]byte(nil), core1.CASAVA...),
					},
					{
						seq:    append([]byte(nil), core2.Seq...),
						qual:   append([]byte(nil), core2.Qual...),
						casava: append([]byte(nil), core2.CASAVA...),
					},
				},
			})
		}
	}

	if len(chunk) > 0 {
		chunks <- chunk
	}
	close(chunks)

	wg.Wait()
	close(results)
	collectWG.Wait()

	defer func() {
		for _, f := range tmpFiles {
			os.Remove(f)
		}
		os.Remove(tmpDir)
	}()

	if readErr != nil {
		return readErr
	}

	select {
	case err := <-errCh:
		return err
	default:
	}

	return mergeFqChunks(tmpFiles, outPrefix, converter, standardize, numMates)
}

// formatConverter picks which fastq.CoreFq writer method reconstructs
// output records: the format-specific one matching fmtType (so output
// round-trips the input format), or ToStandard when standardize is set.
func formatConverter(fmtType string, standardize bool) (func(*fastq.CoreFq, *bufio.Writer), error) {
	if standardize {
		return (*fastq.CoreFq).ToStandard, nil
	}
	switch fmtType {
	case "haplotagging":
		return (*fastq.CoreFq).ToHaplotagging, nil
	case "stlfr":
		return (*fastq.CoreFq).ToStlfr, nil
	case "tellseq":
		return (*fastq.CoreFq).ToTellseq, nil
	case "10x":
		return (*fastq.CoreFq).ToTenX, nil
	default:
		return nil, fmt.Errorf("unsupported linked-read technology %q for format-preserving output; use --standardize", fmtType)
	}
}

// sortAndSpill sorts one chunk by barcode and writes it to a fresh,
// uncompressed temp file in djinn's own binary pair format (see
// writePair/readPair) -- not FASTQ text, since this file is only ever
// read back by mergeFqChunks and re-parsing FASTQ text would mean paying
// the barcode/CASAVA extraction cost a second time for no benefit.
func fqsortAndSpill(chunk []fqPair, tmpDir string) (string, error) {
	slices.SortFunc(chunk, func(a, b fqPair) int { return bytes.Compare(pairKey(&a), pairKey(&b)) })

	f, err := os.CreateTemp(tmpDir, ".fqsort-*.bin")
	if err != nil {
		return "", err
	}
	defer f.Close()

	w := bufio.NewWriterSize(f, 1<<20)
	for i := range chunk {
		if err := writePair(w, &chunk[i]); err != nil {
			return "", err
		}
	}
	if err := w.Flush(); err != nil {
		return "", err
	}
	return f.Name(), nil
}

// ── binary spill format ─────────────────────────────────────────────────────

func writeBytes(w *bufio.Writer, b []byte) error {
	var lenBuf [4]byte
	binary.LittleEndian.PutUint32(lenBuf[:], uint32(len(b)))
	if _, err := w.Write(lenBuf[:]); err != nil {
		return err
	}
	if len(b) > 0 {
		if _, err := w.Write(b); err != nil {
			return err
		}
	}
	return nil
}

func readBytes(r *bufio.Reader) ([]byte, error) {
	var lenBuf [4]byte
	if _, err := io.ReadFull(r, lenBuf[:]); err != nil {
		return nil, err
	}
	n := binary.LittleEndian.Uint32(lenBuf[:])
	if n == 0 {
		return nil, nil
	}
	b := make([]byte, n)
	if _, err := io.ReadFull(r, b); err != nil {
		return nil, err
	}
	return b, nil
}

// flagValid is the one bit packed into writePair/readPair's single
// trailing byte.
const flagValid = 1 << 0

// writePair/readPair's field order must match exactly, and both must
// agree (out of band -- see numMates below) on how many mates each
// record carries. Barcode leads since it's the sort key -- were this
// format ever inspected line-by-line (e.g. hexdump during debugging)
// the key would be the first thing in each record, same spirit as
// putting it first in a hypothetical TSV encoding.
func writePair(w *bufio.Writer, p *fqPair) error {
	for _, b := range [...][]byte{p.barcode, p.id, p.comments} {
		if err := writeBytes(w, b); err != nil {
			return err
		}
	}
	for _, m := range p.mates {
		if err := writeBytes(w, m.seq); err != nil {
			return err
		}
		if err := writeBytes(w, m.qual); err != nil {
			return err
		}
		if err := writeBytes(w, m.casava); err != nil {
			return err
		}
	}
	var flags byte
	if p.valid {
		flags |= flagValid
	}
	return w.WriteByte(flags)
}

// readPair reads one record written by writePair, which carried
// numMates mates (fixed for the whole run -- every spill file from the
// same SortFqByBX call agrees). err is io.EOF (and p is nil) only when
// the stream is exhausted exactly at a record boundary, same contract
// sort.mergeChunks relies on for its BAM readers.
func readPair(r *bufio.Reader, numMates int) (*fqPair, error) {
	barcode, err := readBytes(r)
	if err != nil {
		return nil, err
	}
	id, err := readBytes(r)
	if err != nil {
		return nil, err
	}
	comments, err := readBytes(r)
	if err != nil {
		return nil, err
	}
	mates := make([]mate, numMates)
	for i := range mates {
		seq, err := readBytes(r)
		if err != nil {
			return nil, err
		}
		qual, err := readBytes(r)
		if err != nil {
			return nil, err
		}
		casava, err := readBytes(r)
		if err != nil {
			return nil, err
		}
		mates[i] = mate{seq: seq, qual: qual, casava: casava}
	}
	flags, err := r.ReadByte()
	if err != nil {
		return nil, err
	}
	return &fqPair{
		barcode:  barcode,
		id:       id,
		comments: comments,
		valid:    flags&flagValid != 0,
		mates:    mates,
	}, nil
}

// ── k-way merge ──────────────────────────────────────────────────────────────

// fqMergeItem is one entry in the merge heap: a pair and which temp file
// it came from. Hand-rolled (not container/heap) for the same reason as
// sort.mergeHeap: avoid boxing every pair into an interface{} on every
// heap operation.
type fqMergeItem struct {
	pair *fqPair
	src  int
}

type fqMergeHeap []fqMergeItem

func (h *fqMergeHeap) push(it fqMergeItem) {
	*h = append(*h, it)
	i := len(*h) - 1
	for i > 0 {
		parent := (i - 1) / 2
		if bytes.Compare(pairKey((*h)[i].pair), pairKey((*h)[parent].pair)) >= 0 {
			break
		}
		(*h)[i], (*h)[parent] = (*h)[parent], (*h)[i]
		i = parent
	}
}

func (h *fqMergeHeap) pop() fqMergeItem {
	old := *h
	n := len(old) - 1
	top := old[0]
	old[0] = old[n]
	*h = old[:n]
	h.siftDown(0)
	return top
}

func (h *fqMergeHeap) siftDown(i int) {
	n := len(*h)
	for {
		left := 2*i + 1
		if left >= n {
			return
		}
		smallest := left
		if right := left + 1; right < n && bytes.Compare(pairKey((*h)[right].pair), pairKey((*h)[left].pair)) < 0 {
			smallest = right
		}
		if bytes.Compare(pairKey((*h)[smallest].pair), pairKey((*h)[i].pair)) >= 0 {
			return
		}
		(*h)[i], (*h)[smallest] = (*h)[smallest], (*h)[i]
		i = smallest
	}
}

// mergeFqChunks does a single-threaded k-way merge of the already
// barcode-sorted temp files, writing the result as gzipped FASTQ via
// converter -- one output file (outPrefix+".R1.fq.gz") for single-end,
// two (adding ".R2.fq.gz") for paired, per numMates. standardize controls
// how an absent CASAVA marker is handled -- see writePairAsFastq.
func mergeFqChunks(tmpFiles []string, outPrefix string, converter func(*fastq.CoreFq, *bufio.Writer), standardize bool, numMates int) (err error) {
	readers := make([]*bufio.Reader, len(tmpFiles))
	files := make([]*os.File, len(tmpFiles))
	for i, path := range tmpFiles {
		f, ferr := os.Open(path)
		if ferr != nil {
			return ferr
		}
		files[i] = f
		readers[i] = bufio.NewReaderSize(f, 1<<20)
	}
	defer func() {
		for _, f := range files {
			f.Close()
		}
	}()

	outs := make([]*xopen.Writer, numMates)
	writers := make([]*bufio.Writer, numMates)
	for i := range outs {
		suffix := ".R1.fq.gz"
		if i == 1 {
			suffix = ".R2.fq.gz"
		}
		o, oerr := xopen.Wopen(outPrefix + suffix)
		if oerr != nil {
			return oerr
		}
		outs[i] = o
		writers[i] = o.Writer
	}
	defer func() {
		for _, o := range outs {
			o.Close()
		}
	}()

	h := &fqMergeHeap{}
	for i, r := range readers {
		p, perr := readPair(r, numMates)
		if perr == nil {
			h.push(fqMergeItem{pair: p, src: i})
		} else if perr != io.EOF {
			return perr
		}
	}

	for len(*h) > 0 {
		top := h.pop()
		writePairAsFastq(top.pair, writers, converter, standardize)
		p, perr := readPair(readers[top.src], numMates)
		if perr == nil {
			h.push(fqMergeItem{pair: p, src: top.src})
		} else if perr != io.EOF {
			return perr
		}
	}

	return nil
}

// mateCasavaByte returns the fallback single-byte CASAVA marker for mate
// index i (0 -> '1', 1 -> '2'), used only when a mate genuinely had no
// marker of its own -- see writePairAsFastq.
func mateCasavaByte(i int) []byte {
	if i == 0 {
		return casavaOne
	}
	return casavaTwo
}

// writePairAsFastq re-assembles each mate as a fastq.CoreFq and writes
// it with converter, reproducing either the original format (converter
// is one of ToHaplotagging/ToStlfr/ToTellseq/ToTenX) or the generic
// BX:Z:/VX:i: form every other djinn fastq tool already produces
// (converter is ToStandard, when standardize is set).
//
// A mate's own CASAVA marker, when it has one, is always passed through
// verbatim -- not just its first byte -- since ToStlfr/ToTellseq/ToTenX
// write the whole thing back out for anything longer than one byte (a
// new-style Illumina marker's real filter flag and index sequence,
// which a synthesized placeholder can't reproduce). Only ToStandard
// treats an absent marker as optional (it omits the "/1"-style suffix
// entirely, matching every other CoreFq-producing caller that never
// found one); the format-specific writers always expect *a* mate digit
// to embed when there are two mates -- ToHaplotagging in particular
// indexes CASAVA[0] unconditionally and would panic on a nil/empty
// CASAVA -- so for paired output with standardize false, a mate that
// truly had no marker of its own gets a bare single-byte fallback
// synthesized from which side of the pair it's on. Single-end output
// has no "other side" to distinguish from, so there's nothing to
// synthesize: the marker is only ever written if the source record
// actually had one, same as ToStandard's own rule.
func writePairAsFastq(p *fqPair, writers []*bufio.Writer, converter func(*fastq.CoreFq, *bufio.Writer), standardize bool) {
	forceMarker := !standardize && len(p.mates) == 2
	for i, m := range p.mates {
		casava := m.casava
		if len(casava) == 0 && forceMarker {
			casava = mateCasavaByte(i)
		}
		c := fastq.CoreFq{ID: p.id, Seq: m.seq, Qual: m.qual, CASAVA: casava, Barcode: p.barcode, Valid: p.valid, Comments: p.comments}
		converter(&c, writers[i])
	}
}
