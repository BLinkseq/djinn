package invalid

import (
	"testing"
)

func TestFilterInvalidXam_SplitsValidAndInvalid(t *testing.T) {
	in := tempPath(t, "in.sam")
	writeSamFile(t, in, []samRec{
		{name: "good1", bx: "AAAA", vx: -1},
		{name: "bad1", bx: "NNNN", vx: -1}, // invalid: contains N
		{name: "nobarcode", bx: "", vx: -1},
	})

	invalidOut := tempPath(t, "invalid.sam")

	out, err := captureStdout(t, func() error {
		return FilterInvalidXam(in, invalidOut, true, 1)
	})
	if err != nil {
		t.Fatalf("FilterInvalidXam returned error: %v", err)
	}

	validRecs := parseSamText(t, out)
	validNames := recNames(validRecs)
	if !contains(validNames, "good1") {
		t.Errorf("expected good1 in valid output, got %v", validNames)
	}
	if contains(validNames, "bad1") || contains(validNames, "nobarcode") {
		t.Errorf("valid output should not contain invalid/no-barcode reads, got %v", validNames)
	}

	invalidRecs := readAllSamRecords(t, invalidOut)
	invalidNames := recNames(invalidRecs)
	if !contains(invalidNames, "bad1") {
		t.Errorf("expected bad1 in invalid output, got %v", invalidNames)
	}
	if !contains(invalidNames, "nobarcode") {
		t.Errorf("expected nobarcode in invalid output, got %v", invalidNames)
	}
}

func TestFilterInvalidXam_NoInvalidOutput(t *testing.T) {
	in := tempPath(t, "in.sam")
	writeSamFile(t, in, []samRec{
		{name: "good1", bx: "AAAA", vx: -1},
		{name: "bad1", bx: "NNNN", vx: -1},
	})

	out, err := captureStdout(t, func() error {
		return FilterInvalidXam(in, "", true, 1)
	})
	if err != nil {
		t.Fatalf("FilterInvalidXam returned error: %v", err)
	}
	validNames := recNames(parseSamText(t, out))
	if !contains(validNames, "good1") {
		t.Errorf("expected good1 in valid output, got %v", validNames)
	}
	if contains(validNames, "bad1") {
		t.Errorf("bad1 should not be in valid output, got %v", validNames)
	}
}

func TestFilterInvalidXam_EmptyInput(t *testing.T) {
	in := tempPath(t, "empty.sam")
	writeSamFile(t, in, nil)

	out, err := captureStdout(t, func() error {
		return FilterInvalidXam(in, "", true, 1)
	})
	if err != nil {
		t.Fatalf("FilterInvalidXam returned error on empty input: %v", err)
	}
	if recs := parseSamText(t, out); len(recs) != 0 {
		t.Errorf("expected no records for an empty SAM file, got %d", len(recs))
	}
}
