package ncbi

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// decodeBamRecords parses out as a real BAM stream and returns every
// record. It fails the test outright if out isn't valid BAM -- a bare
// len(out) != 0 check would pass even on a headerless/malformed stream.
func decodeBamRecords(t *testing.T, out []byte) []*sam.Record {
	t.Helper()
	br, err := bam.NewReader(bytes.NewReader(out), 1)
	if err != nil {
		t.Fatalf("output is not valid BAM: %v", err)
	}
	defer br.Close()

	var recs []*sam.Record
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("reading BAM record: %v", err)
		}
		recs = append(recs, rec)
	}
	return recs
}

const ncbiTenXBC = "AAAACCCCGGGGTTTT" // 16nt

func writeNcbiTenXWhitelist(t *testing.T, dir string) string {
	t.Helper()
	path := filepath.Join(dir, "whitelist.txt")
	if err := os.WriteFile(path, []byte(ncbiTenXBC+"\n"), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	return path
}

func writeNcbiTenXFastq(t *testing.T, path, id, seqSuffix string) {
	t.Helper()
	qual := strings.Repeat("I", len(ncbiTenXBC)+len(seqSuffix))
	content := "@" + id + "\n" + ncbiTenXBC + seqSuffix + "\n+\n" + qual + "\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

func writeNcbiPlainFastq(t *testing.T, path, id, seq string) {
	t.Helper()
	qual := strings.Repeat("J", len(seq))
	content := "@" + id + "\n" + seq + "\n+\n" + qual + "\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

// captureStdoutBytes redirects os.Stdout to a pipe for the duration of fn
// and returns everything written to it. fn (NCBIfq10x/NcbiFq) closes
// os.Stdout itself when it succeeds, which is exactly what signals the
// drain goroutine below to stop; os.Stdout is restored to orig
// unconditionally afterward so later tests in this process aren't left
// with a closed stdout.
func captureStdoutBytes(t *testing.T, fn func() error) ([]byte, error) {
	t.Helper()
	orig := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("os.Pipe: %v", err)
	}
	os.Stdout = w

	done := make(chan []byte, 1)
	go func() {
		var buf bytes.Buffer
		io.Copy(&buf, r)
		done <- buf.Bytes()
	}()

	fnErr := fn()
	if fnErr != nil {
		w.Close()
	}
	os.Stdout = orig
	return <-done, fnErr
}

func TestNCBIfq10x_SingleFile(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeNcbiTenXFastq(t, r1, "read1", "ACGT")
	wl := writeNcbiTenXWhitelist(t, dir)

	out, err := captureStdoutBytes(t, func() error {
		return NCBIfq10x([]string{r1}, wl, 1)
	})
	if err != nil {
		t.Fatalf("NCBIfq10x: %v", err)
	}
	recs := decodeBamRecords(t, out)
	if len(recs) != 1 {
		t.Fatalf("got %d records, want 1", len(recs))
	}
	rec := recs[0]
	if rec.Name != "read1" {
		t.Errorf("Name = %q, want %q", rec.Name, "read1")
	}
	if rec.Flags&sam.Unmapped == 0 {
		t.Errorf("Flags = %v, want Unmapped set", rec.Flags)
	}
	if rec.Flags&sam.Paired != 0 {
		t.Errorf("Flags = %v, want Paired unset for single-end input", rec.Flags)
	}
	if got := string(rec.Seq.Expand()); got != "ACGT" {
		t.Errorf("Seq = %q, want %q (barcode must not leak into the stored sequence)", got, "ACGT")
	}
}

func TestNCBIfq10x_PairedFiles(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	writeNcbiTenXFastq(t, r1, "read1/1", "ACGT")
	writeNcbiPlainFastq(t, r2, "read1/2", "TGCA")
	wl := writeNcbiTenXWhitelist(t, dir)

	out, err := captureStdoutBytes(t, func() error {
		return NCBIfq10x([]string{r1, r2}, wl, 1)
	})
	if err != nil {
		t.Fatalf("NCBIfq10x: %v", err)
	}
	recs := decodeBamRecords(t, out)
	if len(recs) != 2 {
		t.Fatalf("got %d records, want 2", len(recs))
	}
	r1rec, r2rec := recs[0], recs[1]
	if r1rec.Flags&sam.Paired == 0 || r1rec.Flags&sam.Read1 == 0 {
		t.Errorf("R1 Flags = %v, want Paired|Read1 set", r1rec.Flags)
	}
	if r2rec.Flags&sam.Paired == 0 || r2rec.Flags&sam.Read2 == 0 {
		t.Errorf("R2 Flags = %v, want Paired|Read2 set", r2rec.Flags)
	}
	if got := string(r1rec.Seq.Expand()); got != "ACGT" {
		t.Errorf("R1 Seq = %q, want %q", got, "ACGT")
	}
	if got := string(r2rec.Seq.Expand()); got != "TGCA" {
		t.Errorf("R2 Seq = %q, want %q", got, "TGCA")
	}
}

func TestNCBIfq10x_MismatchedMateIDsIsError(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	writeNcbiTenXFastq(t, r1, "read1", "ACGT")
	writeNcbiPlainFastq(t, r2, "different", "TGCA")
	wl := writeNcbiTenXWhitelist(t, dir)

	_, err := captureStdoutBytes(t, func() error {
		return NCBIfq10x([]string{r1, r2}, wl, 1)
	})
	if err == nil {
		t.Fatal("expected an error for mismatched mate IDs, got nil")
	}
}

func TestNcbiFq_SingleFile(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeNcbiPlainFastq(t, r1, "read1 BX:Z:A01C02B03D04", "ACGT")

	out, err := captureStdoutBytes(t, func() error {
		return NcbiFq([]string{r1}, "", 1)
	})
	if err != nil {
		t.Fatalf("NcbiFq: %v", err)
	}
	recs := decodeBamRecords(t, out)
	if len(recs) != 1 {
		t.Fatalf("got %d records, want 1", len(recs))
	}
	rec := recs[0]
	if rec.Name != "read1" {
		t.Errorf("Name = %q, want %q", rec.Name, "read1")
	}
	if got := string(rec.Seq.Expand()); got != "ACGT" {
		t.Errorf("Seq = %q, want %q", got, "ACGT")
	}
	var bx *sam.Aux
	for _, a := range rec.AuxFields {
		if a.Tag() == sam.NewTag("BX") {
			bx = &a
			break
		}
	}
	if bx == nil {
		t.Fatalf("expected a BX aux tag, got %v", rec.AuxFields)
	}
	if got := bx.Value(); got != "A01C02B03D04" {
		t.Errorf("BX value = %v, want %q", got, "A01C02B03D04")
	}
}

func TestNcbiFq_PairedFiles(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	writeNcbiPlainFastq(t, r1, "read1 BX:Z:A01C02B03D04", "ACGT")
	writeNcbiPlainFastq(t, r2, "read1 BX:Z:A01C02B03D04", "TGCA")

	out, err := captureStdoutBytes(t, func() error {
		return NcbiFq([]string{r1, r2}, "", 1)
	})
	if err != nil {
		t.Fatalf("NcbiFq: %v", err)
	}
	recs := decodeBamRecords(t, out)
	if len(recs) != 2 {
		t.Fatalf("got %d records, want 2", len(recs))
	}
	if got := string(recs[0].Seq.Expand()); got != "ACGT" {
		t.Errorf("R1 Seq = %q, want %q", got, "ACGT")
	}
	if got := string(recs[1].Seq.Expand()); got != "TGCA" {
		t.Errorf("R2 Seq = %q, want %q", got, "TGCA")
	}
}
