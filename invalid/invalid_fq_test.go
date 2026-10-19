package invalid

import (
	"path/filepath"
	"testing"
)

func TestFilterInvalidFQ_SplitsValidAndInvalid(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")

	reads := []hapRead{
		{id: "good1", bc: "AAAA", vx: 1},
		{id: "bad1", bc: "NNNN", vx: 0},
		{id: "nobarcode", bc: "", vx: -1},
	}
	reads = padHapReads(reads, 100)
	writeHapFastq(t, r1, reads)

	prefix := filepath.Join(dir, "out")
	invPrefix := filepath.Join(dir, "inv")

	if err := FilterInvalidFQ([]string{r1}, prefix, invPrefix); err != nil {
		t.Fatalf("FilterInvalidFQ returned error: %v", err)
	}

	validIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	invalidIDs := readGzFastqIDs(t, invPrefix+".R1.fq.gz")

	if !contains(validIDs, "good1") {
		t.Errorf("expected good1 in valid output, got %v", validIDs)
	}
	if contains(validIDs, "bad1") || contains(validIDs, "nobarcode") {
		t.Errorf("valid output should not contain invalid/no-barcode reads, got %v", validIDs)
	}
	if !contains(invalidIDs, "bad1") {
		t.Errorf("expected bad1 in invalid output, got %v", invalidIDs)
	}
	if !contains(invalidIDs, "nobarcode") {
		t.Errorf("expected nobarcode (no BX tag) to be routed to invalid output, got %v", invalidIDs)
	}
}

func TestFilterInvalidFQ_NoInvalidPrefix(t *testing.T) {
	// When invalidprefix == "", invalid reads should simply be dropped, and
	// no ".invalid" style file should be created.
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	reads := padHapReads([]hapRead{
		{id: "good1", bc: "AAAA", vx: 1},
		{id: "bad1", bc: "NNNN", vx: 0},
	}, 100)
	writeHapFastq(t, r1, reads)

	prefix := filepath.Join(dir, "out")
	if err := FilterInvalidFQ([]string{r1}, prefix, ""); err != nil {
		t.Fatalf("FilterInvalidFQ returned error: %v", err)
	}

	validIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	if !contains(validIDs, "good1") {
		t.Errorf("expected good1 in valid output, got %v", validIDs)
	}
	if contains(validIDs, "bad1") {
		t.Errorf("bad1 should have been dropped, got %v", validIDs)
	}
}

func TestFilterInvalidFQ_PairedFiles(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")

	r1Reads := padHapReads([]hapRead{{id: "p1", bc: "AAAA", vx: 1}}, 100)
	writeHapFastq(t, r1, r1Reads)
	// R2 is short: format detection always keys off fqs[0], so this is fine.
	writeHapFastq(t, r2, []hapRead{{id: "p1", bc: "TTTT", vx: 1}})

	prefix := filepath.Join(dir, "out")
	if err := FilterInvalidFQ([]string{r1, r2}, prefix, ""); err != nil {
		t.Fatalf("FilterInvalidFQ returned error: %v", err)
	}
	if ids := readGzFastqIDs(t, prefix+".R1.fq.gz"); !contains(ids, "p1") {
		t.Errorf("expected p1 in R1 output, got %v", ids)
	}
	if ids := readGzFastqIDs(t, prefix+".R2.fq.gz"); !contains(ids, "p1") {
		t.Errorf("expected p1 in R2 output, got %v", ids)
	}
}

func TestFilterInvalidFQ_MissingFile(t *testing.T) {
	dir := t.TempDir()
	err := FilterInvalidFQ([]string{"/nonexistent/file.fq"}, filepath.Join(dir, "out"), "")
	if err == nil {
		t.Fatal("expected an error for a nonexistent input file")
	}
}

func contains(ss []string, v string) bool {
	for _, s := range ss {
		if s == v {
			return true
		}
	}
	return false
}
