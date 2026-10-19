package ncbi

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const ncbiSamHeader = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"

func writeNcbiSam(t *testing.T, path, content string) {
	t.Helper()
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

func readGzFastq(t *testing.T, path string) string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Open %s: %v", path, err)
	}
	defer f.Close()
	gz, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip.NewReader %s: %v", path, err)
	}
	defer gz.Close()
	var buf bytes.Buffer
	if _, err := io.Copy(&buf, gz); err != nil {
		t.Fatalf("reading gzip content of %s: %v", path, err)
	}
	return buf.String()
}

func TestNCBIRoutesReadsByFlag(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	prefix := filepath.Join(dir, "out")

	// flag 77 = paired|unmapped|mate_unmapped|read1
	// flag 141 = paired|unmapped|mate_unmapped|read2
	samText := ncbiSamHeader +
		"pair1\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"pair1\t141\t*\t0\t0\t*\t*\t0\t0\tTGCATGCATG\tIIIIIIIIII\n"
	writeNcbiSam(t, in, samText)

	if err := NcbiXam(in, prefix, 1); err != nil {
		t.Fatalf("NCBI: %v", err)
	}

	r1Path := prefix + ".R1.fq.gz"
	r2Path := prefix + ".R2.fq.gz"
	if _, err := os.Stat(r1Path); err != nil {
		t.Fatalf("expected R1 output file: %v", err)
	}
	if _, err := os.Stat(r2Path); err != nil {
		t.Fatalf("expected R2 output file: %v", err)
	}

	r1 := readGzFastq(t, r1Path)
	r2 := readGzFastq(t, r2Path)

	if !strings.Contains(r1, "@pair1/1") {
		t.Errorf("R1 output missing expected read header, got:\n%s", r1)
	}
	if !strings.Contains(r1, "ACGTACGTAC") {
		t.Errorf("R1 output missing expected sequence, got:\n%s", r1)
	}
	if !strings.Contains(r2, "@pair1/2") {
		t.Errorf("R2 output missing expected read header, got:\n%s", r2)
	}
	if !strings.Contains(r2, "TGCATGCATG") {
		t.Errorf("R2 output missing expected sequence, got:\n%s", r2)
	}
}

func TestNCBIEmptyInput(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	prefix := filepath.Join(dir, "out")
	writeNcbiSam(t, in, ncbiSamHeader)

	if err := NcbiXam(in, prefix, 1); err != nil {
		t.Fatalf("NCBI: %v", err)
	}

	r1 := readGzFastq(t, prefix+".R1.fq.gz")
	r2 := readGzFastq(t, prefix+".R2.fq.gz")
	if r1 != "" {
		t.Errorf("expected empty R1 output, got %q", r1)
	}
	if r2 != "" {
		t.Errorf("expected empty R2 output, got %q", r2)
	}
}

// TestNCBIUnpairedReadRoutedToR2 documents current (questionable) behavior:
// a record with neither the Read1 nor Read2 flag set (e.g. a genuinely
// single-end/unpaired read) is routed to the R2 output because NCBI() only
// tests `rec.Flags&sam.Read1 != 0` and treats every other record as R2. This
// is flagged in the audit report as a design concern rather than fixed here,
// since the correct behavior for non-paired input is ambiguous (skip? error?
// write to a third file?) and out of scope for a minimal bug fix.
func TestNCBIUnpairedReadRoutedToR2(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	prefix := filepath.Join(dir, "out")

	// flag 4 = unmapped, single-end (no Paired/Read1/Read2 bits at all)
	samText := ncbiSamHeader +
		"single1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	writeNcbiSam(t, in, samText)

	if err := NcbiXam(in, prefix, 1); err != nil {
		t.Fatalf("NCBI: %v", err)
	}

	r1 := readGzFastq(t, prefix+".R1.fq.gz")
	r2 := readGzFastq(t, prefix+".R2.fq.gz")

	if r1 != "" {
		t.Errorf("expected empty R1 output for unpaired read, got %q", r1)
	}
	if !strings.Contains(r2, "@single1/2") {
		t.Errorf("current behavior: unpaired read expected (incorrectly) in R2 output, got R2=%q", r2)
	}
}
