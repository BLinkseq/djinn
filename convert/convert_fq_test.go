package convert

import (
	"bufio"
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const tenXBC = "AAAACCCCGGGGTTTT" // 16nt

func writeTenXWhitelist(t *testing.T, dir string, barcodes ...string) string {
	t.Helper()
	path := filepath.Join(dir, "whitelist.txt")
	content := strings.Join(barcodes, "\n") + "\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	return path
}

func writeTenXFastq(t *testing.T, path, id, seqSuffix string) {
	t.Helper()
	qual := strings.Repeat("I", len(tenXBC)+len(seqSuffix))
	content := "@" + id + "\n" + tenXBC + seqSuffix + "\n+\n" + qual + "\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

func writePlainFastq(t *testing.T, path, id, seq string) {
	t.Helper()
	qual := strings.Repeat("J", len(seq))
	content := "@" + id + "\n" + seq + "\n+\n" + qual + "\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
}

// readGzFastqHeaders reads a gzip FASTQ file and returns each record's
// full header line (minus the leading '@').
func readGzFastqHeaders(t *testing.T, path string) []string {
	t.Helper()
	fh, err := os.Open(path)
	if err != nil {
		t.Fatalf("os.Open(%s): %v", path, err)
	}
	defer fh.Close()
	gz, err := gzip.NewReader(fh)
	if err != nil {
		t.Fatalf("gzip.NewReader: %v", err)
	}
	defer gz.Close()

	var headers []string
	s := bufio.NewScanner(gz)
	for s.Scan() {
		line := s.Text()
		if strings.HasPrefix(line, "@") {
			headers = append(headers, strings.TrimPrefix(line, "@"))
			for range 3 {
				s.Scan()
			}
		}
	}
	return headers
}

func fileExists(path string) bool {
	_, err := os.Stat(path)
	return err == nil
}

// ── ConvertFrom10X ───────────────────────────────────────────────────────────

func TestConvertFrom10X_SingleFile(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeTenXFastq(t, r1, "read1", "ACGT")

	wl := writeTenXWhitelist(t, dir, tenXBC)
	prefix := filepath.Join(dir, "out")

	if err := ConvertFrom10X([]string{r1}, "haplotagging", prefix, wl); err != nil {
		t.Fatalf("ConvertFrom10X: %v", err)
	}

	headers := readGzFastqHeaders(t, prefix+".R1.fq.gz")
	if len(headers) != 1 {
		t.Fatalf("got %d R1 records, want 1", len(headers))
	}
	if !strings.Contains(headers[0], "A01C01B01D01") {
		t.Errorf("R1 header = %q, want it to contain the first-generated haplotagging barcode A01C01B01D01", headers[0])
	}

	if fileExists(prefix + ".R2.fq.gz") {
		t.Errorf("expected no R2 output file for single-end input")
	}
	if !fileExists(prefix + ".bc.map") {
		t.Errorf("expected a .bc.map file to be written")
	}
}

func TestConvertFrom10X_PairedFiles(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	writeTenXFastq(t, r1, "read1", "ACGT")
	writePlainFastq(t, r2, "read1", "TGCA")

	wl := writeTenXWhitelist(t, dir, tenXBC)
	prefix := filepath.Join(dir, "out")

	if err := ConvertFrom10X([]string{r1, r2}, "haplotagging", prefix, wl); err != nil {
		t.Fatalf("ConvertFrom10X: %v", err)
	}

	r1Headers := readGzFastqHeaders(t, prefix+".R1.fq.gz")
	r2Headers := readGzFastqHeaders(t, prefix+".R2.fq.gz")
	if len(r1Headers) != 1 || len(r2Headers) != 1 {
		t.Fatalf("got %d R1 and %d R2 records, want 1 each", len(r1Headers), len(r2Headers))
	}
	// R2 has no barcode of its own; it must carry R1's converted barcode.
	if !strings.Contains(r1Headers[0], "A01C01B01D01") || !strings.Contains(r2Headers[0], "A01C01B01D01") {
		t.Errorf("expected both mates to carry the same converted barcode, got R1=%q R2=%q", r1Headers[0], r2Headers[0])
	}
}

func TestConvertFrom10X_MismatchedMateIDsIsError(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	writeTenXFastq(t, r1, "read1", "ACGT")
	writePlainFastq(t, r2, "different", "TGCA")

	wl := writeTenXWhitelist(t, dir, tenXBC)
	prefix := filepath.Join(dir, "out")

	if err := ConvertFrom10X([]string{r1, r2}, "haplotagging", prefix, wl); err == nil {
		t.Fatal("expected an error for mismatched mate IDs, got nil")
	}
}

// ── StandardizeFrom10X ───────────────────────────────────────────────────────

func TestStandardizeFrom10X_SingleFile(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeTenXFastq(t, r1, "read1", "ACGT")

	wl := writeTenXWhitelist(t, dir, tenXBC)
	prefix := filepath.Join(dir, "out")

	if err := StandardizeFrom10X([]string{r1}, prefix, wl); err != nil {
		t.Fatalf("StandardizeFrom10X: %v", err)
	}

	headers := readGzFastqHeaders(t, prefix+".R1.fq.gz")
	if len(headers) != 1 {
		t.Fatalf("got %d R1 records, want 1", len(headers))
	}
	if !strings.Contains(headers[0], "BX:Z:"+tenXBC) {
		t.Errorf("R1 header = %q, want it to contain BX:Z:%s", headers[0], tenXBC)
	}
	if fileExists(prefix + ".R2.fq.gz") {
		t.Errorf("expected no R2 output file for single-end input")
	}
}

func TestStandardizeFrom10X_PairedFiles(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	writeTenXFastq(t, r1, "read1", "ACGT")
	writePlainFastq(t, r2, "read1", "TGCA")

	wl := writeTenXWhitelist(t, dir, tenXBC)
	prefix := filepath.Join(dir, "out")

	if err := StandardizeFrom10X([]string{r1, r2}, prefix, wl); err != nil {
		t.Fatalf("StandardizeFrom10X: %v", err)
	}

	r1Headers := readGzFastqHeaders(t, prefix+".R1.fq.gz")
	r2Headers := readGzFastqHeaders(t, prefix+".R2.fq.gz")
	if len(r1Headers) != 1 || len(r2Headers) != 1 {
		t.Fatalf("got %d R1 and %d R2 records, want 1 each", len(r1Headers), len(r2Headers))
	}
	if !strings.Contains(r1Headers[0], "BX:Z:"+tenXBC) || !strings.Contains(r2Headers[0], "BX:Z:"+tenXBC) {
		t.Errorf("expected both mates to carry R1's barcode, got R1=%q R2=%q", r1Headers[0], r2Headers[0])
	}
}
