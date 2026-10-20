package sample

import (
	"os"
	"path/filepath"
	"testing"
)

// tenBarcodeReads builds 100 reads across 10 distinct valid barcodes
// (10 reads each), so nBC == 10 exactly.
func tenBarcodeReads() []hapRead {
	var reads []hapRead
	for b := range 10 {
		for i := range 10 {
			reads = append(reads, hapRead{
				id: "r", bc: barcodeName(b), vx: 1,
			})
			_ = i
		}
	}
	return reads
}

func barcodeName(i int) string {
	letters := "ACGTACGTAC"
	return string([]byte{letters[i], letters[i], letters[i], letters[i]}) + string(rune('A'+i))
}

func TestSampleFq_FractionalDownsample(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastq(t, r1, tenBarcodeReads())

	prefix := filepath.Join(dir, "out")
	if err := SampleFq([]string{r1}, prefix, 0.3, 42, false); err != nil {
		t.Fatalf("SampleFq returned error: %v", err)
	}

	bcLines := readLines(t, prefix+".bc")
	if len(bcLines) != 3 {
		t.Fatalf("expected 3 sampled barcodes (round(10*0.3)), got %d: %v", len(bcLines), bcLines)
	}

	// The .bc file must live alongside the prefix (not in the CWD).
	if _, err := os.Stat(prefix + ".bc"); err != nil {
		t.Fatalf(".bc file not found at expected prefixed path: %v", err)
	}
}

func TestSampleFq_AbsoluteCountDownsample(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastq(t, r1, tenBarcodeReads())

	prefix := filepath.Join(dir, "out")
	if err := SampleFq([]string{r1}, prefix, 4, 42, false); err != nil {
		t.Fatalf("SampleFq returned error: %v", err)
	}
	bcLines := readLines(t, prefix+".bc")
	if len(bcLines) != 4 {
		t.Fatalf("expected 4 sampled barcodes, got %d: %v", len(bcLines), bcLines)
	}
}

func TestSampleFq_DeterministicWithSeed(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastq(t, r1, tenBarcodeReads())

	prefixA := filepath.Join(dir, "a")
	prefixB := filepath.Join(dir, "b")

	if err := SampleFq([]string{r1}, prefixA, 5, 7, false); err != nil {
		t.Fatalf("SampleFq (a) returned error: %v", err)
	}
	if err := SampleFq([]string{r1}, prefixB, 5, 7, false); err != nil {
		t.Fatalf("SampleFq (b) returned error: %v", err)
	}

	linesA := readLines(t, prefixA+".bc")
	linesB := readLines(t, prefixB+".bc")
	if len(linesA) != len(linesB) {
		t.Fatalf("sample sizes differ: %d vs %d", len(linesA), len(linesB))
	}
	for i := range linesA {
		if linesA[i] != linesB[i] {
			t.Fatalf("same seed produced different sampled barcode order/content at index %d: %q vs %q", i, linesA[i], linesB[i])
		}
	}
}

func TestSampleFq_TooFewBarcodes(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastq(t, r1, tenBarcodeReads())

	prefix := filepath.Join(dir, "out")
	err := SampleFq([]string{r1}, prefix, 50, 1, false)
	if err == nil {
		t.Fatal("expected an error when downsample count exceeds the number of available barcodes")
	}
}

func TestSampleFq_OutputContainsOnlySampledBarcodes(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastq(t, r1, tenBarcodeReads())

	prefix := filepath.Join(dir, "out")
	if err := SampleFq([]string{r1}, prefix, 3, 1, false); err != nil {
		t.Fatalf("SampleFq returned error: %v", err)
	}
	kept := readLines(t, prefix+".bc")
	keptSet := map[string]bool{}
	for _, k := range kept {
		keptSet[k] = true
	}

	ids := readGzFastqIDs(t, prefix+".R1.fq.gz")
	if len(ids) == 0 {
		t.Fatal("expected some reads in sampled output")
	}
	// All emitted reads carry one of the kept barcodes; since every read
	// shares the id "r", we instead check the record count matches
	// (10 reads per barcode) * (number kept barcodes).
	wantCount := len(kept) * 10
	if len(ids) != wantCount {
		t.Errorf("expected %d reads in sampled output (10 per kept barcode), got %d", wantCount, len(ids))
	}
}

func TestSampleFq_MissingFile(t *testing.T) {
	dir := t.TempDir()
	err := SampleFq([]string{"/nonexistent/file.fq"}, filepath.Join(dir, "out"), 1, 1, false)
	if err == nil {
		t.Fatal("expected an error for a nonexistent input file")
	}
}

func TestSampleFq_FewerThan100Records(t *testing.T) {
	// fastq.CoreFqParser samples up to 100 records to detect the
	// linked-read technology; a file with fewer records should still be
	// handled cleanly rather than failing outright.
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastq(t, r1, []hapRead{{bc: "AAAA", vx: 1}})

	err := SampleFq([]string{r1}, filepath.Join(dir, "out"), 1, 1, false)
	if err != nil {
		t.Fatalf("expected no error for a fastq file with fewer than 100 records, got %v", err)
	}
}
