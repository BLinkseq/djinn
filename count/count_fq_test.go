package count

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestCountFQ_R1Only(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")

	reads := []hapRead{
		{bc: "AAAA", vx: 1},
		{bc: "AAAA", vx: 1},
		{bc: "AAAA", vx: 1},
		{bc: "CCCC", vx: 1},
		{bc: "CCCC", vx: 1},
		{bc: "GGGG", vx: 0}, // invalid
		{bc: "", vx: 1},     // no barcode found -> ignored regardless of tag values
	}
	reads = padHapReads(reads, 100)
	writeHapFastq(t, r1, reads)

	out, err := captureStdout(t, func() error {
		return CountFQ([]string{r1}, false, "")
	})
	if err != nil {
		t.Fatalf("CountFQ returned error: %v", err)
	}
	counts := parseCounts(t, out)

	if counts["AAAA"] != 3 {
		t.Errorf("AAAA count = %d, want 3", counts["AAAA"])
	}
	if counts["CCCC"] != 2 {
		t.Errorf("CCCC count = %d, want 2", counts["CCCC"])
	}
	if _, ok := counts["GGGG"]; ok {
		t.Errorf("GGGG (invalid barcode) should be excluded when invalid=false, got %d", counts["GGGG"])
	}
}

func TestCountFQ_IncludeInvalid(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")

	reads := []hapRead{
		{bc: "GGGG", vx: 0},
		{bc: "GGGG", vx: 0},
	}
	reads = padHapReads(reads, 100)
	writeHapFastq(t, r1, reads)

	out, err := captureStdout(t, func() error {
		return CountFQ([]string{r1}, true, "")
	})
	if err != nil {
		t.Fatalf("CountFQ returned error: %v", err)
	}
	counts := parseCounts(t, out)
	if counts["GGGG"] != 2 {
		t.Errorf("GGGG count = %d, want 2 (invalid=true should include invalid barcodes)", counts["GGGG"])
	}
}

func TestCountFQ_R1R2Union(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")

	r1Reads := padHapReads([]hapRead{
		{bc: "AAAA", vx: 1},
		{bc: "AAAA", vx: 1},
	}, 100)
	writeHapFastq(t, r1, r1Reads)

	// R2 does not need 100 records: only fqs[0] goes through format
	// detection.
	r2Reads := []hapRead{
		{bc: "AAAA", vx: 1}, // already counted via R1 -> must not double count
		{bc: "TTTT", vx: 1}, // new barcode, seen twice in R2
		{bc: "TTTT", vx: 1},
	}
	writeHapFastq(t, r2, r2Reads)

	out, err := captureStdout(t, func() error {
		return CountFQ([]string{r1, r2}, false, "")
	})
	if err != nil {
		t.Fatalf("CountFQ returned error: %v", err)
	}
	counts := parseCounts(t, out)

	if counts["AAAA"] != 2 {
		t.Errorf("AAAA count = %d, want 2 (R2 occurrence of an R1 barcode must not be double counted)", counts["AAAA"])
	}
	if counts["TTTT"] != 2 {
		t.Errorf("TTTT count = %d, want 2", counts["TTTT"])
	}
}

func TestCountFQ_FewerThan100Records(t *testing.T) {
	// fastq.CoreFqParser samples up to 100 records to detect the
	// linked-read technology; a file with fewer records should still be
	// handled cleanly (format detection proceeds with whatever it read)
	// rather than failing outright.
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastq(t, r1, []hapRead{{bc: "AAAA", vx: 1}})

	_, err := captureStdout(t, func() error {
		return CountFQ([]string{r1}, false, "")
	})
	if err != nil {
		t.Fatalf("expected no error for a fastq file with fewer than 100 records, got %v", err)
	}
}

func TestCountFQ_MissingFile(t *testing.T) {
	_, err := captureStdout(t, func() error {
		return CountFQ([]string{"/nonexistent/path/does-not-exist.fq"}, false, "")
	})
	if err == nil {
		t.Fatal("expected an error opening a nonexistent file, got nil")
	}
}

// writeTenXFastq writes a plain-text FASTQ file of 10X-format records: a
// 16bp barcode inline as the sequence prefix, no separate tag.
func writeTenXFastq(t *testing.T, path string, ids []string, barcode string) {
	t.Helper()
	var sb strings.Builder
	for _, id := range ids {
		seq := barcode + "ACGTACGTAC"
		sb.WriteString("@" + id + "\n" + seq + "\n+\n" + strings.Repeat("I", len(seq)) + "\n")
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

func writeTenXWhitelist(t *testing.T, path string, barcodes ...string) {
	t.Helper()
	if err := os.WriteFile(path, []byte(strings.Join(barcodes, "\n")+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}
}

func TestCount10XFQ(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	wl := filepath.Join(dir, "whitelist.txt")

	writeTenXWhitelist(t, wl, "AAAACCCCGGGGTTTT")
	// two reads with the known barcode, one with an unrecognized prefix
	writeTenXFastq(t, r1, []string{"read1", "read2"}, "AAAACCCCGGGGTTTT")
	writeTenXFastq(t, filepath.Join(dir, "r1b.fq"), []string{"read3"}, "ZZZZZZZZZZZZZZZZ")
	// append the unrecognized-barcode read to r1 so both are in one file
	data, err := os.ReadFile(r1)
	if err != nil {
		t.Fatal(err)
	}
	extra, err := os.ReadFile(filepath.Join(dir, "r1b.fq"))
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(r1, append(data, extra...), 0o644); err != nil {
		t.Fatal(err)
	}

	out, err := captureStdout(t, func() error {
		return CountFQ([]string{r1}, false, wl)
	})
	if err != nil {
		t.Fatalf("CountFQ (10X): %v", err)
	}
	counts := parseCounts(t, out)
	if counts["AAAACCCCGGGGTTTT"] != 2 {
		t.Errorf("AAAACCCCGGGGTTTT count = %d, want 2", counts["AAAACCCCGGGGTTTT"])
	}
	if len(counts) != 1 {
		t.Errorf("expected only the whitelisted barcode to be counted, got %v", counts)
	}
}

func TestCount10XFQ_EmbeddedWhitelistSentinel(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeTenXFastq(t, r1, []string{"read1"}, "AAAACCCCGGGGTTTT")

	// "10x"/"10X" means "use the embedded official whitelist" (fastq.Load10X),
	// not a literal file named "10x" -- this used to be unreachable because
	// the CLI's own Args validator ran filecheck("10x") before this point,
	// but CountFQ/Count10XFQ themselves have always honored the sentinel.
	_, err := captureStdout(t, func() error {
		return CountFQ([]string{r1}, false, "10x")
	})
	if err != nil {
		t.Fatalf("CountFQ with embedded whitelist sentinel: %v", err)
	}
}
