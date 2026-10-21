package unlinked

import (
	"os"
	"path/filepath"
	"slices"
	"strings"
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

	if err := RemoveUnlinkedFQ([]string{r1}, prefix, true, true, ""); err != nil {
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

	if _, err := os.Stat(prefix + ".count"); err != nil {
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
	if err := RemoveUnlinkedFQ([]string{r1}, prefix, false, false, ""); err != nil {
		t.Fatalf("RemoveUnlinkedFQ returned error: %v", err)
	}
	mainIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	if !slices.Contains(mainIDs, "a1") || slices.Contains(mainIDs, "c1") {
		t.Errorf("unexpected main output contents: %v", mainIDs)
	}
}

func TestRemoveUnlinkedFQ_MissingFile(t *testing.T) {
	dir := t.TempDir()
	err := RemoveUnlinkedFQ([]string{"/nonexistent/file.fq"}, filepath.Join(dir, "out"), false, false, "")
	if err == nil {
		t.Fatal("expected an error for a nonexistent input file")
	}
}

func writeTenXR1Unlinked(t *testing.T, path string, reads [][2]string) {
	t.Helper()
	var sb strings.Builder
	for _, r := range reads {
		id, barcode := r[0], r[1]
		seq := barcode + "ACGTACGTAC"
		sb.WriteString("@" + id + "\n" + seq + "\n+\n" + strings.Repeat("I", len(seq)) + "\n")
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

func writeTenXR2Unlinked(t *testing.T, path string, ids []string) {
	t.Helper()
	var sb strings.Builder
	for _, id := range ids {
		sb.WriteString("@" + id + "\nGGGGCCCCAAAATTTT\n+\nIIIIIIIIIIIIIIII\n")
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

// Regression test: RemoveUnlinkedFQ accepted a bcmap parameter but never
// used it anywhere in the function body, so --bc-list/-b was silently a
// no-op for rm-unlinked -- a genuine 10X input would always fall through
// to auto-detection and fail with "unable to determine linked-read
// technology", since 10X has no self-describing marker to detect.
func TestRemoveUnlinked10XFQ_SplitsSingletonsAndNot(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	wl := filepath.Join(dir, "whitelist.txt")
	if err := os.WriteFile(wl, []byte("AAAACCCCGGGGTTTT\nTTTTGGGGCCCCAAAA\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	writeTenXR1Unlinked(t, r1, [][2]string{
		{"a1", "AAAACCCCGGGGTTTT"},
		{"a2", "AAAACCCCGGGGTTTT"}, // same barcode as a1 -> linked
		{"c1", "TTTTGGGGCCCCAAAA"}, // singleton
	})
	writeTenXR2Unlinked(t, r2, []string{"a1", "a2", "c1"})

	prefix := filepath.Join(dir, "out")
	if err := RemoveUnlinkedFQ([]string{r1, r2}, prefix, true, true, wl); err != nil {
		t.Fatalf("RemoveUnlinkedFQ (10X): %v", err)
	}

	mainR1 := readGzFastqIDs(t, prefix+".R1.fq.gz")
	mainR2 := readGzFastqIDs(t, prefix+".R2.fq.gz")
	singleR1 := readGzFastqIDs(t, prefix+".unlinked.R1.fq.gz")
	singleR2 := readGzFastqIDs(t, prefix+".unlinked.R2.fq.gz")

	hasPrefix := func(ids []string, id string) bool {
		return slices.ContainsFunc(ids, func(s string) bool { return strings.HasPrefix(s, id) })
	}

	if !hasPrefix(mainR1, "a1") || !hasPrefix(mainR1, "a2") {
		t.Errorf("expected a1,a2 in main R1 output, got %v", mainR1)
	}
	if !hasPrefix(mainR2, "a1") || !hasPrefix(mainR2, "a2") {
		t.Errorf("expected a1,a2 in main R2 output, got %v", mainR2)
	}
	if hasPrefix(mainR1, "c1") || hasPrefix(mainR2, "c1") {
		t.Errorf("singleton c1 should not be in main output, got R1=%v R2=%v", mainR1, mainR2)
	}
	if !hasPrefix(singleR1, "c1") || !hasPrefix(singleR2, "c1") {
		t.Errorf("expected c1 in singleton output, got R1=%v R2=%v", singleR1, singleR2)
	}
	if _, err := os.Stat(prefix + ".count"); err != nil {
		t.Errorf("expected barcode count file to be created: %v", err)
	}
}

func TestRemoveUnlinked10XFQ_SingleEnd(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	wl := filepath.Join(dir, "whitelist.txt")
	if err := os.WriteFile(wl, []byte("AAAACCCCGGGGTTTT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	writeTenXR1Unlinked(t, r1, [][2]string{
		{"a1", "AAAACCCCGGGGTTTT"},
		{"a2", "AAAACCCCGGGGTTTT"},
	})

	prefix := filepath.Join(dir, "out")
	if err := RemoveUnlinkedFQ([]string{r1}, prefix, false, false, wl); err != nil {
		t.Fatalf("RemoveUnlinkedFQ (10X single-end): %v", err)
	}
	if _, err := os.Stat(prefix + ".R2.fq.gz"); err == nil {
		t.Error("single-end 10X input must not produce an R2 output file")
	}
	mainR1 := readGzFastqIDs(t, prefix+".R1.fq.gz")
	if len(mainR1) != 2 {
		t.Errorf("got %d records in main output, want 2: %v", len(mainR1), mainR1)
	}
}
