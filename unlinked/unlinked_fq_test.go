package unlinked

import (
	"os"
	"path/filepath"
	"slices"
	"testing"
)

func TestGetFqCount_Basic(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")

	reads := []hapRead{
		{id: "a1", bc: "AAAA", vx: 1},
		{id: "a2", bc: "AAAA", vx: 1},
		{id: "a3", bc: "AAAA", vx: 1},
		{id: "c1", bc: "CCCC", vx: 1}, // singleton
		{id: "g1", bc: "GGGG", vx: 0}, // invalid, must not be counted
		{id: "g2", bc: "GGGG", vx: 0},
		{id: "n1", bc: "", vx: -1}, // no barcode
	}
	reads = padHapReads(reads, 100)
	writeHapFastq(t, r1, reads)

	counts, err := getFqCount([]string{r1})
	if err != nil {
		t.Fatalf("getFqCount returned error: %v", err)
	}
	if counts["AAAA"] != 3 {
		t.Errorf("AAAA count = %d, want 3", counts["AAAA"])
	}
	if counts["CCCC"] != 1 {
		t.Errorf("CCCC count = %d, want 1", counts["CCCC"])
	}
	if _, ok := counts["GGGG"]; ok {
		t.Errorf("invalid barcode GGGG must not be counted, got %d", counts["GGGG"])
	}
}

func TestGetFqCount_R1R2Union(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")

	r1Reads := padHapReads([]hapRead{
		{bc: "AAAA", vx: 1},
	}, 100)
	writeHapFastq(t, r1, r1Reads)

	r2Reads := []hapRead{
		{bc: "AAAA", vx: 1}, // must not double count (already in R1 set)
		{bc: "TTTT", vx: 1},
		{bc: "TTTT", vx: 1},
	}
	writeHapFastq(t, r2, r2Reads)

	counts, err := getFqCount([]string{r1, r2})
	if err != nil {
		t.Fatalf("getFqCount returned error: %v", err)
	}
	if counts["AAAA"] != 1 {
		t.Errorf("AAAA count = %d, want 1 (R2 occurrence must not double count an R1 barcode)", counts["AAAA"])
	}
	if counts["TTTT"] != 2 {
		t.Errorf("TTTT count = %d, want 2", counts["TTTT"])
	}
}

func TestGetFqCount_MissingFile(t *testing.T) {
	_, err := getFqCount([]string{"/nonexistent/file.fq"})
	if err == nil {
		t.Fatal("expected an error for a nonexistent input file")
	}
}

func TestRemoveUnlinkedFQ_SplitsSingletonsAndNot(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")

	reads := []hapRead{
		{id: "a1", bc: "AAAA", vx: 1},
		{id: "a2", bc: "AAAA", vx: 1},
		{id: "c1", bc: "CCCC", vx: 1}, // singleton
		{id: "g1", bc: "GGGG", vx: 0}, // invalid, dropped entirely
	}
	reads = padHapReads(reads, 100)
	writeHapFastq(t, r1, reads)

	prefix := filepath.Join(dir, "out")
	bcCountFile := filepath.Join(dir, "counts.tsv")

	if err := RemoveUnlinkedFQ([]string{r1}, prefix, true, bcCountFile); err != nil {
		t.Fatalf("RemoveUnlinkedFQ returned error: %v", err)
	}

	mainIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	singleIDs := readGzFastqIDs(t, prefix+".unlinked.R1.fq.gz")

	if !slices.Contains(mainIDs, "a1") || !slices.Contains(mainIDs, "a2") {
		t.Errorf("expected a1,a2 (non-singleton barcode) in main output, got %v", mainIDs)
	}
	if slices.Contains(mainIDs, "c1") {
		t.Errorf("singleton read c1 should not be in main output, got %v", mainIDs)
	}
	if !slices.Contains(singleIDs, "c1") {
		t.Errorf("expected c1 (singleton barcode) in singleton output, got %v", singleIDs)
	}
	if slices.Contains(mainIDs, "g1") || slices.Contains(singleIDs, "g1") {
		t.Errorf("invalid-barcode read g1 should be dropped entirely, got main=%v single=%v", mainIDs, singleIDs)
	}

	if _, err := os.Stat(bcCountFile); err != nil {
		t.Errorf("expected barcode count file to be created: %v", err)
	}
}

func TestRemoveUnlinkedFQ_NoSingletonOutput(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	reads := padHapReads([]hapRead{
		{id: "a1", bc: "AAAA", vx: 1},
		{id: "a2", bc: "AAAA", vx: 1},
		{id: "c1", bc: "CCCC", vx: 1},
	}, 100)
	writeHapFastq(t, r1, reads)

	prefix := filepath.Join(dir, "out")
	if err := RemoveUnlinkedFQ([]string{r1}, prefix, false, ""); err != nil {
		t.Fatalf("RemoveUnlinkedFQ returned error: %v", err)
	}
	mainIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	if !slices.Contains(mainIDs, "a1") || slices.Contains(mainIDs, "c1") {
		t.Errorf("unexpected main output contents: %v", mainIDs)
	}
}

func TestRemoveUnlinkedFQ_MissingFile(t *testing.T) {
	dir := t.TempDir()
	err := RemoveUnlinkedFQ([]string{"/nonexistent/file.fq"}, filepath.Join(dir, "out"), false, "")
	if err == nil {
		t.Fatal("expected an error for a nonexistent input file")
	}
}
