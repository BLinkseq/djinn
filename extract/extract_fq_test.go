package extract

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeFastq(t *testing.T, path string, records []string) {
	t.Helper()
	var sb strings.Builder
	for _, r := range records {
		sb.WriteString(r)
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

func haplotagRecord(name, desc, seqStr string) string {
	return "@" + name + " " + desc + "\n" + seqStr + "\n+\n" + strings.Repeat("I", len(seqStr)) + "\n"
}

func TestExtractFQHaplotagging(t *testing.T) {
	dir := t.TempDir()
	fq := filepath.Join(dir, "reads.fq")
	writeFastq(t, fq, []string{
		haplotagRecord("read1", "BX:Z:A01C01B01D01", "ACGTACGTAC"),
		haplotagRecord("read2", "BX:Z:A01C01B01D01", "ACGTACGTAC"), // duplicate barcode
		haplotagRecord("read3", "BX:Z:A00C01B01D01", "ACGTACGTAC"), // invalid barcode
		haplotagRecord("read4", "no barcode here", "ACGTACGTAC"),   // no barcode
	})

	var err error
	out := captureStdout(t, func() {
		err = ExtractFQ([]string{fq}, false)
	})
	if err != nil {
		t.Fatalf("ExtractFQ: %v", err)
	}

	lines := splitNonEmpty(string(out))
	if len(lines) != 1 {
		t.Fatalf("got %d barcodes, want 1: %v", len(lines), lines)
	}
	if lines[0] != "A01C01B01D01" {
		t.Errorf("barcode = %q, want %q", lines[0], "A01C01B01D01")
	}
}

func TestExtractFQIncludeInvalid(t *testing.T) {
	dir := t.TempDir()
	fq := filepath.Join(dir, "reads.fq")
	writeFastq(t, fq, []string{
		haplotagRecord("read1", "BX:Z:A01C01B01D01", "ACGTACGTAC"),
		haplotagRecord("read2", "BX:Z:A00C01B01D01", "ACGTACGTAC"), // invalid barcode
	})

	var err error
	out := captureStdout(t, func() {
		err = ExtractFQ([]string{fq}, true)
	})
	if err != nil {
		t.Fatalf("ExtractFQ: %v", err)
	}
	lines := splitNonEmpty(string(out))
	if len(lines) != 2 {
		t.Fatalf("got %d barcodes, want 2 (invalid included): %v", len(lines), lines)
	}
}

func TestExtractFQMultipleFiles(t *testing.T) {
	dir := t.TempDir()
	fq1 := filepath.Join(dir, "r1.fq")
	fq2 := filepath.Join(dir, "r2.fq")
	writeFastq(t, fq1, []string{
		haplotagRecord("read1", "BX:Z:A01C01B01D01", "ACGTACGTAC"),
	})
	writeFastq(t, fq2, []string{
		haplotagRecord("read2", "BX:Z:A02C02B02D02", "ACGTACGTAC"),
		haplotagRecord("read3", "BX:Z:A01C01B01D01", "ACGTACGTAC"), // duplicate of fq1's barcode
	})

	var err error
	out := captureStdout(t, func() {
		err = ExtractFQ([]string{fq1, fq2}, false)
	})
	if err != nil {
		t.Fatalf("ExtractFQ: %v", err)
	}
	lines := splitNonEmpty(string(out))
	if len(lines) != 2 {
		t.Fatalf("got %d barcodes across files, want 2 (deduped across files): %v", len(lines), lines)
	}
}

func TestExtractFQBadFirstFile(t *testing.T) {
	err := ExtractFQ([]string{"/nonexistent/path/reads.fq"}, false)
	if err == nil {
		t.Fatal("expected error for nonexistent input file, got nil")
	}
}
