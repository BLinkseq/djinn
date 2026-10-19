package stats

import (
	"fmt"
	"os"
	"slices"
	"strings"
	"testing"
)

// ── golden outputs from the Python original ───────────────────────────────────
//
// testdata/fq/<platform>.tsv is what harpy's bx_stats_fq.py prints for testdata/fq/<platform>.fq
// (see testdata/regen.sh). The output order is fixed, so it is compared exactly.

func TestStatsFQ_GoldenPython(t *testing.T) {
	for _, platform := range []string{Haplotagging, Stlfr, Tellseq} {
		want, err := os.ReadFile("testdata/fq/" + platform + ".tsv")
		if err != nil {
			t.Fatal(err)
		}
		// the platform named explicitly, in any capitalization, and detected from the reads
		for _, given := range []string{platform, strings.ToUpper(platform), ""} {
			got, err := captureStdout(t, func() error { return StatsFQ("testdata/fq/"+platform+".fq", given) })
			if err != nil {
				t.Fatalf("%s (platform %q): %v", platform, given, err)
			}
			if got != string(want) {
				t.Errorf("%s (platform %q): output differs from the Python version\ngot:\n%s\nwant:\n%s", platform, given, got, want)
			}
		}
	}
}

// ── helpers ───────────────────────────────────────────────────────────────────

// fq builds one FASTQ record. header is everything after the "@" on the first line.
func fq(header string) string { return "@" + header + "\nACGTACGTAC\n+\nIIIIIIIIII\n" }

func writeFastq(t *testing.T, records ...string) string {
	t.Helper()
	path := tempPath(t, "in.fq")
	if err := os.WriteFile(path, []byte(strings.Join(records, "")), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

func statsFQOf(t *testing.T, path, platform string) string {
	t.Helper()
	out, err := captureStdout(t, func() error { return StatsFQ(path, platform) })
	if err != nil {
		t.Fatalf("StatsFQ: %v", err)
	}
	return out
}

// ── which positions of a barcode are bad ──────────────────────────────────────

func TestSchemeMark(t *testing.T) {
	for _, c := range []struct {
		name   string
		scheme scheme
		bc     string
		want   []int
	}{
		{"haplotagging, all good", haplotaggingScheme, "A01C02B03D04", []int{0, 0, 0, 0}},
		{"haplotagging, A and B", haplotaggingScheme, "A00C02B00D04", []int{1, 0, 1, 0}},
		{"haplotagging, all four", haplotaggingScheme, "A00C00B00D00", []int{1, 1, 1, 1}},
		{"haplotagging, 00 inside a larger number is not a failed segment", haplotaggingScheme, "A10C02B03D04", []int{0, 0, 0, 0}},
		{"haplotagging, lowercase letters", haplotaggingScheme, "a00c02b03d00", []int{1, 0, 0, 1}},
		{"haplotagging, empty", haplotaggingScheme, "", []int{0, 0, 0, 0}},
		{"stlfr, all good", stlfrScheme, "1_2_3", []int{0, 0, 0}},
		{"stlfr, first and last", stlfrScheme, "0_5_0", []int{1, 0, 1}},
		{"stlfr, middle", stlfrScheme, "10_0_30", []int{0, 1, 0}},
		{"stlfr, 00 or 10 are not zero", stlfrScheme, "00_10_0", []int{0, 0, 1}},
		{"stlfr, extra parts are ignored, not a crash", stlfrScheme, "0_1_2_0", []int{1, 0, 0}},
		{"tellseq, no N", tellseqScheme, strings.Repeat("A", 18), make([]int, 18)},
		{"tellseq, N at 1 and 18", tellseqScheme, "N" + strings.Repeat("C", 16) + "N", append(append([]int{1}, make([]int, 16)...), 1)},
		{"tellseq, a longer barcode lengthens the table", tellseqScheme, strings.Repeat("A", 19) + "N", append(make([]int, 19), 1)},
	} {
		bad := make([]int, c.scheme.slots)
		got := c.scheme.mark([]byte(c.bc), bad)
		if !slices.Equal(got, c.want) {
			t.Errorf("%s: mark(%q) = %v, want %v", c.name, c.bc, got, c.want)
		}
	}
}

func TestSchemeLabels(t *testing.T) {
	var got []string
	for i := range 4 {
		got = append(got, haplotaggingScheme.label(i))
	}
	for i := range 3 {
		got = append(got, stlfrScheme.label(i))
	}
	got = append(got, tellseqScheme.label(0), tellseqScheme.label(17))
	if want := "A C B D 1 2 3 1 18"; strings.Join(got, " ") != want {
		t.Errorf("labels = %v, want %s", got, want)
	}
}

// ── behaviours specific to the djinn version ──────────────────────────────────

// A barcode is valid if the VX tag says so, otherwise if the barcode itself is; the positions that
// read 00 are counted from the barcode text either way.
func TestStatsFQ_ValidityFollowsVXThenBarcode(t *testing.T) {
	path := writeFastq(t,
		fq("r1 BX:Z:A01C02B03D04"),                // no VX, good barcode -> valid
		fq("r2 BX:Z:A00C02B03D04"),                // no VX, A00 -> invalid
		fq("r3 BX:Z:A00C02B03D04\tVX:i:1"),        // VX says valid, so valid, but the A00 is still counted
		fq("r4 BX:Z:A01C02B03D04\tVX:i:0"),        // VX says invalid
		fq("r5 1:N:0:ATGGAGA\tBX:Z:A01C02B00D04"), // CASAVA field ahead of the tag; B00, no VX -> invalid
	)
	got := statsFQOf(t, path, Haplotagging)
	want := "Reads\t5\nBarcodes\t5\nBarcodes_Valid\t2\nBarcodes_Invalid\t3\nA\t2\nC\t0\nB\t1\nD\t0\n" // valid: r1 and r3
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// A read with no comment is simply a read without a barcode (the Python version crashes on it).
func TestStatsFQ_ReadsWithoutBarcode(t *testing.T) {
	path := writeFastq(t,
		fq("r1 BX:Z:A01C02B03D04"),
		fq("r2"),               // no comment at all
		fq("r3 1:N:0:ATGGAGA"), // a comment, but no barcode in it
		fq("r4 BX:Z:"),         // an empty barcode does not count
	)
	got := statsFQOf(t, path, Haplotagging)
	want := "Reads\t4\nBarcodes\t1\nBarcodes_Valid\t1\nBarcodes_Invalid\t0\nA\t0\nC\t0\nB\t0\nD\t0\n"
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// A trailing "-1" (a GEM well number or similar) is not part of the barcode's positions.
func TestStatsFQ_HyphenSuffixIgnored(t *testing.T) {
	path := writeFastq(t,
		fq("r1 BX:Z:A00C02B03D04-1"),
		fq("r2 BX:Z:A01C02B03D04-1"),
	)
	got := statsFQOf(t, path, Haplotagging)
	want := "Reads\t2\nBarcodes\t2\nBarcodes_Valid\t1\nBarcodes_Invalid\t1\nA\t1\nC\t0\nB\t0\nD\t0\n"
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// The suffix must not disturb the positions of the other technologies either: it would hide a final
// stLFR 0 ("0-1" is not "0") and would stretch a TELLseq barcode past 18 bases.
func TestStatsFQ_HyphenSuffixOtherPlatforms(t *testing.T) {
	got := statsFQOf(t, writeFastq(t, fq("a#0_5_0-1"), fq("b#1_2_3-1")), Stlfr)
	if want := "Reads\t2\nBarcodes\t2\nBarcodes_Valid\t1\nBarcodes_Invalid\t1\n1\t1\n2\t0\n3\t1\n"; got != want {
		t.Errorf("stlfr: got:\n%s\nwant:\n%s", got, want)
	}

	got = statsFQOf(t, writeFastq(t, fq("x:1:2:NCGTACGTACGTACGTAC-1")), Tellseq)
	if rows := strings.Count(got, "\n"); rows != 4+18 || !strings.Contains(got, "\n1\t1\n") {
		t.Errorf("tellseq: want 18 position rows with N at 1, got:\n%s", got)
	}
}

func TestStatsFQ_StlfrAndTellseqByName(t *testing.T) {
	got := statsFQOf(t, writeFastq(t, fq("a#1_2_3/1"), fq("b#0_2_0/1"), fq("c/1")), "")
	want := "Reads\t3\nBarcodes\t2\nBarcodes_Valid\t1\nBarcodes_Invalid\t1\n1\t1\n2\t0\n3\t1\n"
	if got != want {
		t.Errorf("stlfr: got:\n%s\nwant:\n%s", got, want)
	}

	got = statsFQOf(t, writeFastq(t, fq("x:1:2:ACGTACGTACGTACGTAC"), fq("y:1:2:NCGTACGTACGTACGTAN"), fq("z:1:2")), "")
	counts := strings.Split(strings.TrimSuffix(got, "\n"), "\n")
	if counts[0] != "Reads\t3" || counts[1] != "Barcodes\t2" || counts[2] != "Barcodes_Valid\t1" ||
		counts[4] != "1\t1" || counts[21] != "18\t1" || len(counts) != 22 {
		t.Errorf("tellseq: unexpected output:\n%s", got)
	}
}

func TestStatsFQ_Errors(t *testing.T) {
	noBarcodes := writeFastq(t, fq("r1 1:N:0:ATGGAGA"), fq("r2 1:N:0:ATGGAGA"))

	if _, err := captureStdout(t, func() error { return StatsFQ(noBarcodes, "") }); err == nil ||
		!strings.Contains(err.Error(), "unable to determine") {
		t.Errorf("no detectable technology: err = %v", err)
	}
	// naming the platform skips the detection: a file with no barcodes is a table of zeros, as in the Python version
	got := statsFQOf(t, noBarcodes, Haplotagging)
	if want := "Reads\t2\nBarcodes\t0\nBarcodes_Valid\t0\nBarcodes_Invalid\t0\nA\t0\nC\t0\nB\t0\nD\t0\n"; got != want {
		t.Errorf("no barcodes with a named platform: got:\n%s\nwant:\n%s", got, want)
	}

	if _, err := captureStdout(t, func() error { return StatsFQ(noBarcodes, "10x") }); err == nil ||
		!strings.Contains(err.Error(), "unknown platform") {
		t.Errorf("unknown platform: err = %v", err)
	}
	if err := StatsFQ(fmt.Sprintf("%s/missing.fq", t.TempDir()), Haplotagging); err == nil {
		t.Error("a missing file should be an error")
	}
}
