package fastq

import (
	"bufio"
	"errors"
	"fmt"
	"os"

	"github.com/biogo/hts/sam"
	"github.com/klauspost/pgzip"
)

type FastqWriter struct {
	f   *os.File
	gz  *pgzip.Writer
	buf *bufio.Writer
	dir string
}

func NewFastqWriter(path string, readnum string, level, threads int) (*FastqWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	gz, err := pgzip.NewWriterLevel(f, level)
	if err != nil {
		f.Close()
		return nil, err
	}
	gz.SetConcurrency(2<<20, threads)
	return &FastqWriter{
		f:   f,
		gz:  gz,
		buf: bufio.NewWriterSize(gz, 2<<20),
		dir: readnum,
	}, nil
}

func (fw *FastqWriter) Close() error {
	var errs []error

	if err := fw.buf.Flush(); err != nil {
		errs = append(errs, fmt.Errorf("flushing buffer: %w", err))
	}
	if err := fw.gz.Close(); err != nil {
		errs = append(errs, fmt.Errorf("closing gzip writer: %w", err))
	}
	if err := fw.f.Close(); err != nil {
		errs = append(errs, fmt.Errorf("closing file: %w", err))
	}

	return errors.Join(errs...)
}

// writeRecord writes a single FASTQ record directly into the bufio buffer.
// qual is raw Phred scores (0-40) as stored in sam.Record — +33 applied inline.
// bufio only flushes to pgzip when its 1MB buffer is full, so no intermediate
// per-record copy is needed.
func (fw *FastqWriter) WriteRecord(name string, auxFields sam.AuxFields, seq []byte, qual []uint8) error {
	w := errWriter{buf: fw.buf}
	w.writeByte('@')
	w.writeString(name)
	w.writeString(fw.dir)
	for _, aux := range auxFields {
		w.writeByte('\t')
		// aux is the raw binary-encoded field (2-byte tag + type + value,
		// no separators); the FASTQ header is text, so it needs the
		// "TAG:TYPE:VALUE" rendering, not the packed bytes.
		w.writeString(aux.String())
	}
	w.writeByte('\n')
	w.write(seq)
	w.writeString("\n+\n")
	if w.err != nil {
		return w.err
	}
	for _, q := range qual {
		if err := fw.buf.WriteByte(q + 33); err != nil {
			return err
		}
	}
	return fw.buf.WriteByte('\n')
}

// errWriter wraps bufio.Writer and stops writing after the first error,
// so call sites can batch writes and check once at the end.
type errWriter struct {
	buf *bufio.Writer
	err error
}

func (w *errWriter) write(b []byte) {
	if w.err == nil {
		_, w.err = w.buf.Write(b)
	}
}

func (w *errWriter) writeString(s string) {
	if w.err == nil {
		_, w.err = w.buf.WriteString(s)
	}
}

func (w *errWriter) writeByte(b byte) {
	if w.err == nil {
		w.err = w.buf.WriteByte(b)
	}
}

func (w *errWriter) writeInt(n int) {
	if w.err == nil {
		w.err = writeInt(w.buf, n)
	}
}

// writeInt writes a non-negative integer to b without allocating.
func writeInt(b *bufio.Writer, n int) error {
	if n == 0 {
		return b.WriteByte('0')
	}
	var tmp [10]byte
	i := len(tmp)
	for n > 0 {
		i--
		tmp[i] = byte('0' + n%10)
		n /= 10
	}
	_, err := b.Write(tmp[i:])
	return err
}
