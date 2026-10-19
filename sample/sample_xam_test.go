package sample

import (
	"os"
	"path/filepath"
	"testing"
)

func tenBarcodeSamRecs() []samRec {
	var recs []samRec
	n := 0
	for b := range 10 {
		for range 10 {
			recs = append(recs, samRec{name: barcodeSamName(n), bx: barcodeName(b), vx: -1})
			n++
		}
	}
	return recs
}

func barcodeSamName(i int) string {
	return "r" + string(rune('0'+i%10)) + string(rune('a'+i/10))
}

func TestSampleXam_FractionalDownsample(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	writeSamFile(t, in, tenBarcodeSamRecs())

	// The .bc file is written to the CWD (basename-derived); isolate it in
	// a scratch directory rather than polluting the real working directory.
	t.Chdir(dir)

	out, err := captureStdout(t, func() error {
		return SampleXam(in, 0.3, 42, 1, false, true)
	})
	if err != nil {
		t.Fatalf("SampleXam returned error: %v", err)
	}

	bcLines := readLines(t, filepath.Base(in)+".bc")
	if len(bcLines) != 3 {
		t.Fatalf("expected 3 sampled barcodes, got %d: %v", len(bcLines), bcLines)
	}

	recs := parseSamText(t, out)
	wantCount := len(bcLines) * 10
	if len(recs) != wantCount {
		t.Errorf("expected %d records in sampled output (10 per kept barcode), got %d", wantCount, len(recs))
	}
}

func TestSampleXam_AbsoluteCountDownsample(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	writeSamFile(t, in, tenBarcodeSamRecs())
	t.Chdir(dir)

	out, err := captureStdout(t, func() error {
		return SampleXam(in, 4, 42, 1, false, true)
	})
	if err != nil {
		t.Fatalf("SampleXam returned error: %v", err)
	}
	bcLines := readLines(t, filepath.Base(in)+".bc")
	if len(bcLines) != 4 {
		t.Fatalf("expected 4 sampled barcodes, got %d: %v", len(bcLines), bcLines)
	}
	recs := parseSamText(t, out)
	if len(recs) != 40 {
		t.Errorf("expected 40 records in output, got %d", len(recs))
	}
}

func TestSampleXam_DeterministicWithSeed(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	writeSamFile(t, in, tenBarcodeSamRecs())
	t.Chdir(dir)

	outA, errA := captureStdout(t, func() error {
		return SampleXam(in, 5, 7, 1, false, true)
	})
	if errA != nil {
		t.Fatalf("SampleXam (a) returned error: %v", errA)
	}
	bcA := append([]string(nil), readLines(t, filepath.Base(in)+".bc")...)

	if err := os.Remove(filepath.Base(in) + ".bc"); err != nil {
		t.Fatal(err)
	}

	outB, errB := captureStdout(t, func() error {
		return SampleXam(in, 5, 7, 1, false, true)
	})
	if errB != nil {
		t.Fatalf("SampleXam (b) returned error: %v", errB)
	}
	bcB := readLines(t, filepath.Base(in)+".bc")

	if len(bcA) != len(bcB) {
		t.Fatalf("sample sizes differ: %d vs %d", len(bcA), len(bcB))
	}
	for i := range bcA {
		if bcA[i] != bcB[i] {
			t.Fatalf("same seed produced different sampled barcodes at index %d: %q vs %q", i, bcA[i], bcB[i])
		}
	}
	if outA != outB {
		t.Errorf("same seed should reproduce identical output stream")
	}
}

func TestSampleXam_TooFewBarcodes(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	writeSamFile(t, in, tenBarcodeSamRecs())
	t.Chdir(dir)

	err := SampleXam(in, 50, 1, 1, false, true)
	if err == nil {
		t.Fatal("expected an error when downsample count exceeds the number of available barcodes")
	}
}

func TestSampleXam_EmptyInput(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "empty.sam")
	writeSamFile(t, in, nil)
	t.Chdir(dir)

	err := SampleXam(in, 1, 1, 1, false, true)
	if err == nil {
		t.Fatal("expected an error when there are 0 barcodes and downsample requests at least 1")
	}
}
