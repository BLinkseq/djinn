package deconvolve

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/biogo/hts/sam"
)

// capture redirects the given *os.File (os.Stdout or os.Stderr) while fn runs and returns what was written to it.
func capture(t *testing.T, target **os.File, fn func() error) (string, error) {
	t.Helper()
	orig := *target
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	*target = w
	outCh := make(chan string, 1)
	go func() {
		var buf bytes.Buffer
		io.Copy(&buf, r)
		outCh <- buf.String()
	}()

	fnErr := fn()

	w.Close()
	*target = orig
	return <-outCh, fnErr
}

const samHeader = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000000\n@SQ\tSN:c2\tLN:1000000\n"

// samLine is one record with a 10M alignment (or the given cigar/flag).
func samLine(name string, flag int, contig string, pos int, cigar string, tags string) string {
	seqLen := 10 // an alignment-less record (cigar "*") still needs some sequence
	if cigar != "*" {
		c, err := sam.ParseCigar([]byte(cigar))
		if err != nil {
			panic(err)
		}
		_, seqLen = c.Lengths() // the sequence must be as long as the CIGAR says
	}
	line := fmt.Sprintf("%s\t%d\t%s\t%d\t60\t%s\t*\t0\t0\t%s\t%s", name, flag, contig, pos, cigar,
		strings.Repeat("A", seqLen), strings.Repeat("I", seqLen))
	if tags != "" {
		line += "\t" + tags
	}
	return line + "\n"
}

// read is an ordinary 10M read.
func read(name, contig string, pos int, tags string) string {
	return samLine(name, 0, contig, pos, "10M", tags)
}

func writeSam(t *testing.T, records ...string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "in.sam")
	if err := os.WriteFile(path, []byte(samHeader+strings.Join(records, "")), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// run deconvolves path with the given distance (writing SAM) and returns the records' (name, BX) pairs,
// the whole stderr summary, and any error. The records come back in file order.
func run(t *testing.T, path string, distance int) (barcodes []string, summary string, err error) {
	t.Helper()
	var out string
	summary, err = capture(t, &os.Stderr, func() error {
		var runErr error
		out, runErr = capture(t, &os.Stdout, func() error { return DeconvolveXam(path, distance, true, 2) })
		return runErr
	})
	for _, line := range strings.Split(out, "\n") {
		if line == "" || strings.HasPrefix(line, "@") {
			continue
		}
		f := strings.Split(line, "\t")
		bx := "-"
		for _, tag := range f[11:] {
			if strings.HasPrefix(tag, "BX:Z:") {
				bx = strings.TrimPrefix(tag, "BX:Z:")
			}
		}
		barcodes = append(barcodes, f[0]+"="+bx)
	}
	return barcodes, summary, err
}
