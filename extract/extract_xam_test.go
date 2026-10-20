package extract

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// captureStdout redirects os.Stdout for the duration of fn and returns
// everything written to it. Extract always writes to "-" (stdout).
func captureStdout(t *testing.T, fn func()) []byte {
	t.Helper()
	orig := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("os.Pipe: %v", err)
	}
	os.Stdout = w

	done := make(chan []byte)
	go func() {
		var buf bytes.Buffer
		io.Copy(&buf, r)
		done <- buf.Bytes()
	}()

	fn()

	w.Close()
	os.Stdout = orig
	return <-done
}

const extractSamHeader = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"

func writeSamFile(t *testing.T, path, content string) {
	t.Helper()
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

func TestExtractValidOnly(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	samText := extractSamHeader +
		"read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A01C01B01D01\n" +
		"read2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A01C01B01D01\n" + // duplicate barcode
		"read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A00C01B01D01\n" + // invalid barcode
		"read4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n" // no barcode
	writeSamFile(t, in, samText)

	var err error
	out := captureStdout(t, func() {
		err = Extract(in, false, 1)
	})
	if err != nil {
		t.Fatalf("Extract: %v", err)
	}

	lines := splitNonEmpty(string(out))
	if len(lines) != 1 {
		t.Fatalf("got %d barcodes, want 1 (dedup, invalid excluded): %v", len(lines), lines)
	}
	if lines[0] != "A01C01B01D01" {
		t.Errorf("barcode = %q, want %q", lines[0], "A01C01B01D01")
	}
}

func TestExtractIncludeInvalid(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	samText := extractSamHeader +
		"read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A01C01B01D01\n" +
		"read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBX:Z:A00C01B01D01\n" // invalid barcode
	writeSamFile(t, in, samText)

	var err error
	out := captureStdout(t, func() {
		err = Extract(in, true, 1)
	})
	if err != nil {
		t.Fatalf("Extract: %v", err)
	}

	lines := splitNonEmpty(string(out))
	if len(lines) != 2 {
		t.Fatalf("got %d barcodes, want 2 (invalid included): %v", len(lines), lines)
	}
}

func TestExtractNoBarcodes(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	samText := extractSamHeader +
		"read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
	writeSamFile(t, in, samText)

	var err error
	out := captureStdout(t, func() {
		err = Extract(in, true, 1)
	})
	if err != nil {
		t.Fatalf("Extract: %v", err)
	}
	if strings.TrimSpace(string(out)) != "" {
		t.Errorf("expected no output, got %q", out)
	}
}

func splitNonEmpty(s string) []string {
	var out []string
	for _, l := range strings.Split(s, "\n") {
		if l != "" {
			out = append(out, l)
		}
	}
	return out
}
