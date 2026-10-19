package count

import (
	"testing"
)

func TestCountXam_Basic(t *testing.T) {
	// "Invalid" is expressed here via barcode content (contains 'N')
	// rather than an explicit VX tag, to keep this test independent of
	// xam.FindBarcode's VX-vs-inferred-validity precedence rules.
	sam := tempPath(t, "in.sam")
	writeSamFile(t, sam, []samRec{
		{name: "r1", bx: "AAAA", vx: -1},
		{name: "r2", bx: "AAAA", vx: -1},
		{name: "r3", bx: "CCCC", vx: -1},
		{name: "r4", bx: "GGGN", vx: -1}, // invalid: contains N
		{name: "r5", bx: "", vx: -1},     // no barcode
	})

	out, err := captureStdout(t, func() error {
		return CountXam(sam, false, 1)
	})
	if err != nil {
		t.Fatalf("CountXam returned error: %v", err)
	}
	counts := parseCounts(t, out)

	if counts["AAAA"] != 2 {
		t.Errorf("AAAA count = %d, want 2", counts["AAAA"])
	}
	if counts["CCCC"] != 1 {
		t.Errorf("CCCC count = %d, want 1", counts["CCCC"])
	}
	if _, ok := counts["GGGN"]; ok {
		t.Errorf("GGGN (invalid) should be excluded when invalid=false")
	}
}

func TestCountXam_IncludeInvalid(t *testing.T) {
	sam := tempPath(t, "in.sam")
	writeSamFile(t, sam, []samRec{
		{name: "r1", bx: "GGGN", vx: -1},
		{name: "r2", bx: "GGGN", vx: -1},
	})

	out, err := captureStdout(t, func() error {
		return CountXam(sam, true, 1)
	})
	if err != nil {
		t.Fatalf("CountXam returned error: %v", err)
	}
	counts := parseCounts(t, out)
	if counts["GGGN"] != 2 {
		t.Errorf("GGGN count = %d, want 2", counts["GGGN"])
	}
}

func TestCountXam_EmptyFile(t *testing.T) {
	sam := tempPath(t, "empty.sam")
	writeSamFile(t, sam, nil)

	out, err := captureStdout(t, func() error {
		return CountXam(sam, false, 1)
	})
	if err != nil {
		t.Fatalf("CountXam returned error on empty input: %v", err)
	}
	if out != "" {
		t.Errorf("expected no output for an empty SAM file, got %q", out)
	}
}
