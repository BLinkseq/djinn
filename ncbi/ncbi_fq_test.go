package ncbi

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

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
	if len(out) == 0 {
		t.Error("expected non-empty BAM output for single-end 10X input")
	}
}

func TestNCBIfq10x_PairedFiles(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	writeNcbiTenXFastq(t, r1, "read1", "ACGT")
	writeNcbiPlainFastq(t, r2, "read1", "TGCA")
	wl := writeNcbiTenXWhitelist(t, dir)

	out, err := captureStdoutBytes(t, func() error {
		return NCBIfq10x([]string{r1, r2}, wl, 1)
	})
	if err != nil {
		t.Fatalf("NCBIfq10x: %v", err)
	}
	if len(out) == 0 {
		t.Error("expected non-empty BAM output for paired 10X input")
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
	if len(out) == 0 {
		t.Error("expected non-empty BAM output for single-end haplotagging input")
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
	if len(out) == 0 {
		t.Error("expected non-empty BAM output for paired haplotagging input")
	}
}
