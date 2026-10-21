package invalid

import (
	"os"
	"path/filepath"
	"slices"
	"strings"
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
	invPrefix := true

	if err := FilterInvalidFQ([]string{r1}, prefix, invPrefix, ""); err != nil {
		t.Fatalf("FilterInvalidFQ returned error: %v", err)
	}

	validIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	invalidIDs := readGzFastqIDs(t, prefix+".invalid.R1.fq.gz")

	if !slices.Contains(validIDs, "good1") {
		t.Errorf("expected good1 in valid output, got %v", validIDs)
	}
	if slices.Contains(validIDs, "bad1") || slices.Contains(validIDs, "nobarcode") {
		t.Errorf("valid output should not contain invalid/no-barcode reads, got %v", validIDs)
	}
	if !slices.Contains(invalidIDs, "bad1") {
		t.Errorf("expected bad1 in invalid output, got %v", invalidIDs)
	}
	if !slices.Contains(invalidIDs, "nobarcode") {
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
	if err := FilterInvalidFQ([]string{r1}, prefix, false, ""); err != nil {
		t.Fatalf("FilterInvalidFQ returned error: %v", err)
	}

	validIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	if !slices.Contains(validIDs, "good1") {
		t.Errorf("expected good1 in valid output, got %v", validIDs)
	}
	if slices.Contains(validIDs, "bad1") {
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
	if err := FilterInvalidFQ([]string{r1, r2}, prefix, false, ""); err != nil {
		t.Fatalf("FilterInvalidFQ returned error: %v", err)
	}
	if ids := readGzFastqIDs(t, prefix+".R1.fq.gz"); !slices.Contains(ids, "p1") {
		t.Errorf("expected p1 in R1 output, got %v", ids)
	}
	if ids := readGzFastqIDs(t, prefix+".R2.fq.gz"); !slices.Contains(ids, "p1") {
		t.Errorf("expected p1 in R2 output, got %v", ids)
	}
}

func TestFilterInvalidFQ_MissingFile(t *testing.T) {
	dir := t.TempDir()
	err := FilterInvalidFQ([]string{"/nonexistent/file.fq"}, filepath.Join(dir, "out"), false, "")
	if err == nil {
		t.Fatal("expected an error for a nonexistent input file")
	}
}

func writeTenXR1(t *testing.T, path string, reads [][2]string) {
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

func writeTenXR2(t *testing.T, path string, ids []string) {
	t.Helper()
	var sb strings.Builder
	for _, id := range ids {
		sb.WriteString("@" + id + "\nGGGGCCCCAAAATTTT\n+\nIIIIIIIIIIIIIIII\n")
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

// Regression test: FilterInvalid10XFQ's paired branch used to never write
// a valid R2 record to its own output file at all (only the keepInvalid
// path wrote R2, for invalid records) -- prefix.R2.fq.gz came out empty
// for every properly-paired, valid-barcode input.
func TestFilterInvalid10XFQ_PairedWritesValidR2(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	wl := filepath.Join(dir, "whitelist.txt")
	if err := os.WriteFile(wl, []byte("AAAACCCCGGGGTTTT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	writeTenXR1(t, r1, [][2]string{{"good1", "AAAACCCCGGGGTTTT"}})
	writeTenXR2(t, r2, []string{"good1"})

	prefix := filepath.Join(dir, "out")
	if err := FilterInvalidFQ([]string{r1, r2}, prefix, false, wl); err != nil {
		t.Fatalf("FilterInvalidFQ (10X paired): %v", err)
	}

	r1IDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	r2IDs := readGzFastqIDs(t, prefix+".R2.fq.gz")
	if !slices.Contains(r1IDs, "good1:AAAACCCCGGGGTTTT") {
		t.Errorf("expected good1 in R1 output, got %v", r1IDs)
	}
	if len(r2IDs) == 0 {
		t.Fatal("R2 output is empty -- valid R2 records must be written, not just invalid ones")
	}
	if !slices.Contains(r2IDs, "good1:AAAACCCCGGGGTTTT") {
		t.Errorf("expected good1 in R2 output, got %v", r2IDs)
	}
}

// Regression test: FilterInvalid10XFQ used to have no separate invalid-output
// file at all for 10X input, unlike every other format -- "kept" invalid
// records were written into the same file as valid ones instead of their
// own prefix.invalid.R*.fq.gz.
func TestFilterInvalid10XFQ_SeparatesInvalidOutput(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	wl := filepath.Join(dir, "whitelist.txt")
	if err := os.WriteFile(wl, []byte("AAAACCCCGGGGTTTT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	writeTenXR1(t, r1, [][2]string{
		{"good1", "AAAACCCCGGGGTTTT"},
		{"bad1", "ZZZZZZZZZZZZZZZZ"},
	})
	writeTenXR2(t, r2, []string{"good1", "bad1"})

	prefix := filepath.Join(dir, "out")
	if err := FilterInvalidFQ([]string{r1, r2}, prefix, true, wl); err != nil {
		t.Fatalf("FilterInvalidFQ (10X paired, keepInvalid): %v", err)
	}

	validR1 := readGzFastqIDs(t, prefix+".R1.fq.gz")
	validR2 := readGzFastqIDs(t, prefix+".R2.fq.gz")
	invalidR1 := readGzFastqIDs(t, prefix+".invalid.R1.fq.gz")
	invalidR2 := readGzFastqIDs(t, prefix+".invalid.R2.fq.gz")

	for _, ids := range [][]string{validR1, validR2} {
		if slices.ContainsFunc(ids, func(s string) bool { return strings.HasPrefix(s, "bad1") }) {
			t.Errorf("invalid record bad1 leaked into the main valid output: %v", ids)
		}
	}
	for _, ids := range [][]string{invalidR1, invalidR2} {
		if !slices.ContainsFunc(ids, func(s string) bool { return strings.HasPrefix(s, "bad1") }) {
			t.Errorf("expected bad1 in the separate invalid output, got %v", ids)
		}
		if slices.ContainsFunc(ids, func(s string) bool { return strings.HasPrefix(s, "good1") }) {
			t.Errorf("valid record good1 leaked into the invalid output: %v", ids)
		}
	}
}

func TestFilterInvalid10XFQ_SingleEnd(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	wl := filepath.Join(dir, "whitelist.txt")
	if err := os.WriteFile(wl, []byte("AAAACCCCGGGGTTTT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	writeTenXR1(t, r1, [][2]string{
		{"good1", "AAAACCCCGGGGTTTT"},
		{"bad1", "ZZZZZZZZZZZZZZZZ"},
	})

	prefix := filepath.Join(dir, "out")
	if err := FilterInvalidFQ([]string{r1}, prefix, true, wl); err != nil {
		t.Fatalf("FilterInvalidFQ (10X single-end, keepInvalid): %v", err)
	}
	if _, err := os.Stat(prefix + ".R2.fq.gz"); err == nil {
		t.Error("single-end 10X input must not produce an R2 output file")
	}
	validIDs := readGzFastqIDs(t, prefix+".R1.fq.gz")
	invalidIDs := readGzFastqIDs(t, prefix+".invalid.R1.fq.gz")
	if !slices.ContainsFunc(validIDs, func(s string) bool { return strings.HasPrefix(s, "good1") }) {
		t.Errorf("expected good1 in valid output, got %v", validIDs)
	}
	if !slices.ContainsFunc(invalidIDs, func(s string) bool { return strings.HasPrefix(s, "bad1") }) {
		t.Errorf("expected bad1 in invalid output, got %v", invalidIDs)
	}
}
