package unlinked

import (
	"os"
	"strings"
	"testing"
)

func TestGetCount_Basic(t *testing.T) {
	in := tempPath(t, "in.sam")
	writeSamFile(t, in, []samRec{
		{name: "a1", bx: "AAAA", vx: -1},
		{name: "a2", bx: "AAAA", vx: -1},
		{name: "c1", bx: "CCCC", vx: -1}, // singleton
		{name: "g1", bx: "GGGN", vx: -1}, // invalid: contains N, must not be counted
		{name: "n1", bx: "", vx: -1},
	})

	counts := getCount(in, 1)
	if counts["AAAA"] != 2 {
		t.Errorf("AAAA count = %d, want 2", counts["AAAA"])
	}
	if counts["CCCC"] != 1 {
		t.Errorf("CCCC count = %d, want 1", counts["CCCC"])
	}
	if _, ok := counts["GGGN"]; ok {
		t.Errorf("invalid barcode must not be counted, got %d", counts["GGGN"])
	}
}

func TestRemoveUnlinkedXam_SplitsSingletonsAndNot(t *testing.T) {
	in := tempPath(t, "in.sam")
	writeSamFile(t, in, []samRec{
		{name: "a1", bx: "AAAA", vx: -1},
		{name: "a2", bx: "AAAA", vx: -1},
		{name: "c1", bx: "CCCC", vx: -1}, // singleton
		{name: "g1", bx: "GGGN", vx: -1}, // invalid, dropped entirely
	})

	singletonOut := tempPath(t, "singletons.sam")
	bcCountFile := tempPath(t, "counts.tsv")

	out, err := captureStdout(t, func() error {
		return RemoveUnlinkedXam(in, singletonOut, bcCountFile, true, 1)
	})
	if err != nil {
		t.Fatalf("RemoveUnlinkedXam returned error: %v", err)
	}

	mainNames := recNames(parseSamText(t, out))
	if !contains(mainNames, "a1") || !contains(mainNames, "a2") {
		t.Errorf("expected a1,a2 (non-singleton) in main output, got %v", mainNames)
	}
	if contains(mainNames, "c1") || contains(mainNames, "g1") {
		t.Errorf("main output should not contain singleton or invalid reads, got %v", mainNames)
	}

	singleNames := recNames(readAllSamRecords(t, singletonOut))
	if !contains(singleNames, "c1") {
		t.Errorf("expected c1 (singleton barcode) in singleton output, got %v", singleNames)
	}
	if contains(singleNames, "g1") {
		t.Errorf("invalid-barcode read g1 should be dropped entirely, not routed to singleton output, got %v", singleNames)
	}

	if _, err := os.Stat(bcCountFile); err != nil {
		t.Errorf("expected barcode count file to be created: %v", err)
	}
}

func TestRemoveUnlinkedXam_NoSingletonOutput(t *testing.T) {
	in := tempPath(t, "in.sam")
	writeSamFile(t, in, []samRec{
		{name: "a1", bx: "AAAA", vx: -1},
		{name: "a2", bx: "AAAA", vx: -1},
		{name: "c1", bx: "CCCC", vx: -1},
	})

	out, err := captureStdout(t, func() error {
		return RemoveUnlinkedXam(in, "", "", true, 1)
	})
	if err != nil {
		t.Fatalf("RemoveUnlinkedXam returned error: %v", err)
	}
	mainNames := recNames(parseSamText(t, out))
	if !contains(mainNames, "a1") || contains(mainNames, "c1") {
		t.Errorf("unexpected main output contents: %v", mainNames)
	}
}

func TestRemoveUnlinkedXam_StdinSpool(t *testing.T) {
	// Exercises the infile == "-" path, which spools stdin to a temp file
	// before reading it twice (once for counts, once for filtering).
	var samText strings.Builder
	samText.WriteString("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n")
	samText.WriteString("a1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tBX:Z:AAAA\n")
	samText.WriteString("a2\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tBX:Z:AAAA\n")
	samText.WriteString("c1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tBX:Z:CCCC\n")

	origStdin := os.Stdin
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	os.Stdin = r
	go func() {
		w.WriteString(samText.String())
		w.Close()
	}()
	t.Cleanup(func() { os.Stdin = origStdin })

	out, err := captureStdout(t, func() error {
		return RemoveUnlinkedXam("-", "", "", true, 1)
	})
	if err != nil {
		t.Fatalf("RemoveUnlinkedXam(\"-\", ...) returned error: %v", err)
	}
	mainNames := recNames(parseSamText(t, out))
	if !contains(mainNames, "a1") || !contains(mainNames, "a2") {
		t.Errorf("expected a1,a2 in main output, got %v", mainNames)
	}
	if contains(mainNames, "c1") {
		t.Errorf("singleton c1 should not be in main output, got %v", mainNames)
	}
}

func TestRemoveUnlinkedXam_EmptyInput(t *testing.T) {
	in := tempPath(t, "empty.sam")
	writeSamFile(t, in, nil)

	out, err := captureStdout(t, func() error {
		return RemoveUnlinkedXam(in, "", "", true, 1)
	})
	if err != nil {
		t.Fatalf("RemoveUnlinkedXam returned error on empty input: %v", err)
	}
	if recs := parseSamText(t, out); len(recs) != 0 {
		t.Errorf("expected no records for an empty SAM file, got %d", len(recs))
	}
}
