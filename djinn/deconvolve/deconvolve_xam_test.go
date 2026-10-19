package deconvolve

import (
	"errors"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"djinn/xam"
)

const aa = "BX:Z:AA\tVX:i:1"

func TestDeconvolve_SuffixesLaterMolecules(t *testing.T) {
	path := writeSam(t,
		read("a", "c1", 100, aa),
		read("b", "c1", 120, aa), // gap 10, same molecule
		read("x", "c1", 200, "BX:Z:BB\tVX:i:1"),
		read("c", "c1", 5000, aa),  // gap far past 100: a new molecule, AA-1
		read("d", "c1", 5100, aa),  // 90 after c's end, still AA-1
		read("e", "c1", 10000, aa), // another jump: AA-2
	)
	got, summary, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"a=AA", "b=AA", "x=BB", "c=AA-1", "d=AA-1", "e=AA-2"}
	if !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
	// barcodes AA (split twice) and BB (never): 2 molecules before, 4 after
	wantSummary := "Original_Molecules\t2\nNew_Molecules\t4\nDeconvolutions\tBarcodes\n0\t1\n2\t1\n"
	if summary != wantSummary {
		t.Errorf("summary:\n%q\nwant:\n%q", summary, wantSummary)
	}
}

// A gap of exactly the distance keeps the molecule together; one base more splits it.
func TestDeconvolve_GapEqualToDistanceDoesNotSplit(t *testing.T) {
	path := writeSam(t,
		read("a", "c1", 100, aa), // ends at 109
		read("b", "c1", 210, aa), // starts at 209: gap exactly 100
		read("c", "c1", 321, aa), // b ends at 219, this starts at 320: gap 101
	)
	got, _, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"a=AA", "b=AA", "c=AA-1"}; !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
}

// The gap is measured from where the previous read ended, not from the furthest end of the molecule so far.
func TestDeconvolve_GapMeasuredFromPreviousRead(t *testing.T) {
	path := writeSam(t,
		samLine("a", 0, "c1", 100, "100M", aa), // ends at 199
		samLine("b", 0, "c1", 150, "10M", aa),  // ends at 159
		samLine("c", 0, "c1", 300, "10M", aa),  // 140 past b, but only 100 past a
	)
	got, _, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"a=AA", "b=AA", "c=AA-1"}; !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
}

func TestDeconvolve_NumberingRestartsOnEachChromosome(t *testing.T) {
	path := writeSam(t,
		read("a", "c1", 100, aa), read("b", "c1", 5000, aa),
		read("c", "c2", 100, aa), read("d", "c2", 5000, aa),
	)
	got, summary, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"a=AA", "b=AA-1", "c=AA", "d=AA-1"}; !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
	// AA on c1 and AA on c2 are two molecules before, each split once
	if want := "Original_Molecules\t2\nNew_Molecules\t4\nDeconvolutions\tBarcodes\n1\t2\n"; summary != want {
		t.Errorf("summary:\n%q\nwant:\n%q", summary, want)
	}
}

// Only valid-barcode reads are touched; the rest of every record, and every other record, is passed through.
func TestDeconvolve_OtherRecordsUnchanged(t *testing.T) {
	path := writeSam(t,
		read("a", "c1", 100, aa),
		read("bad", "c1", 200, "BX:Z:AA\tVX:i:0"),   // invalid by VX
		read("none", "c1", 300, ""),                 // no barcode
		samLine("unmapped", 4, "c1", 400, "*", aa),  // not aligned anywhere
		samLine("unm2", 4, "c1", 3000, "10M", aa),   // flagged unmapped although it has a CIGAR: skipped, it must not split AA
		read("b", "c1", 5000, aa),                   // AA-1
		read("bad2", "c1", 6000, "BX:Z:AA\tVX:i:0"), // still invalid and untouched, even after a split
		read("c", "c1", 10000, "BX:Z:AA-1\tVX:i:1"), // a barcode that already has a hyphen suffix
	)
	got, _, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"a=AA", "bad=AA", "none=-", "unmapped=AA", "unm2=AA", "b=AA-1", "bad2=AA", "c=AA-1"}
	// the last read is the first of its own barcode "AA-1", so it is not split from itself
	if !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
}

// An existing hyphen suffix (like "-1" on a 10X barcode) is kept and the new number added after it.
func TestDeconvolve_ExistingHyphenSuffixKept(t *testing.T) {
	const gem = "BX:Z:ACGT-1\tVX:i:1"
	path := writeSam(t, read("a", "c1", 100, gem), read("b", "c1", 5000, gem))
	got, _, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"a=ACGT-1", "b=ACGT-1-1"}; !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
}

// Supplementary, duplicate and secondary records take the label of their barcode's molecule, but never
// move the barcode along. Once a molecule is closed (the scan has gone more than the distance past its
// last read), a record for that barcode belongs to the next molecule, so it is not counted towards the closed one.
func TestDeconvolve_SupplementaryDuplicateAndSecondary(t *testing.T) {
	path := writeSam(t,
		read("a", "c1", 100, aa),                                   // AA starts; ends at 109
		samLine("sup1", 2048, "c1", 105, "10M", aa),                // while AA is open: AA
		samLine("dup", 1024, "c1", 106, "10M", aa),                 // AA
		samLine("sec", 256, "c1", 107, "10M", aa),                  // AA
		read("b", "c1", 500, "BX:Z:BB\tVX:i:1"),                    // the scan moves 400 past AA's last read: AA's molecule is closed
		samLine("sup2", 2048, "c1", 510, "10M", aa),                // so this belongs to the next molecule: AA-1
		samLine("sup3", 2048, "c1", 520, "10M", "BX:Z:CC\tVX:i:1"), // a barcode not seen yet: first molecule, unchanged
		// AA's next molecule starts here (far from AA's last read, which ended at 109): AA-1. It would be AA-2 if the
		// supplementary above (ending at 519) had been counted as a read, since this is more than 100 past that too.
		read("f", "c1", 640, aa),
	)
	got, _, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"a=AA", "sup1=AA", "dup=AA", "sec=AA", "b=BB", "sup2=AA-1", "sup3=CC", "f=AA-1"}
	if !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
}

// A barcode that comes from the read name has no BX tag; the numbered barcode is written as a BX tag.
func TestDeconvolve_ReadNameBarcodeGetsBXTag(t *testing.T) {
	path := writeSam(t,
		read("r1:ACGTACGTACGTACGTAC", "c1", 100, ""),
		read("r2:ACGTACGTACGTACGTAC", "c1", 5000, ""),
	)
	got, _, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"r1:ACGTACGTACGTACGTAC=-", "r2:ACGTACGTACGTACGTAC=ACGTACGTACGTACGTAC-1"}; !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
}

func TestDeconvolve_NotSorted(t *testing.T) {
	path := writeSam(t, read("a", "c1", 100, aa), read("b", "c1", 50, aa))
	got, summary, err := run(t, path, 100)
	if !errors.Is(err, ErrNotSorted) {
		t.Errorf("err = %v, want ErrNotSorted", err)
	}
	if want := []string{"a=AA"}; !slices.Equal(got, want) {
		t.Errorf("records written before the error = %v, want %v", got, want)
	}
	if summary != "" {
		t.Errorf("no summary should be printed after an error, got %q", summary)
	}
}

func TestDeconvolve_DistanceMustBePositive(t *testing.T) {
	path := writeSam(t, read("a", "c1", 100, aa))
	for _, d := range []int{0, -5} {
		if _, _, err := run(t, path, d); err == nil {
			t.Errorf("distance %d should be an error", d)
		}
	}
}

func TestDeconvolve_EmptyFile(t *testing.T) {
	got, summary, err := run(t, writeSam(t), 100)
	if err != nil || len(got) != 0 {
		t.Fatalf("empty input: records %v, err %v", got, err)
	}
	if want := "Original_Molecules\t0\nNew_Molecules\t0\nDeconvolutions\tBarcodes\n"; summary != want {
		t.Errorf("summary = %q, want %q", summary, want)
	}
}

// The header gains a PG line saying what was done, and BAM output works too.
func TestDeconvolve_HeaderAndBAMOutput(t *testing.T) {
	path := writeSam(t, read("a", "c1", 100, aa), read("b", "c1", 5000, aa))

	out, err := capture(t, &os.Stdout, func() error {
		_, err := capture(t, &os.Stderr, func() error { return DeconvolveXam(path, 100, true, 2) })
		return err
	})
	if err != nil || !strings.Contains(out, "@PG\tID:djinn") || !strings.Contains(out, "CL:djinn sam deconvolve -d 100") {
		t.Errorf("missing @PG line (err %v):\n%s", err, out)
	}

	bamOut, err := capture(t, &os.Stdout, func() error {
		_, err := capture(t, &os.Stderr, func() error { return DeconvolveXam(path, 100, false, 2) })
		return err
	})
	if err != nil {
		t.Fatal(err)
	}
	bamPath := filepath.Join(t.TempDir(), "out.bam")
	if err := os.WriteFile(bamPath, []byte(bamOut), 0o644); err != nil {
		t.Fatal(err)
	}
	recs, _ := xam.NewXamReaderChan(bamPath, 10, xam.IoBuf, 1)
	var bx []string
	for rec := range recs {
		v, _ := xam.GetStringTag(rec, "BX")
		bx = append(bx, v)
	}
	if want := []string{"AA", "AA-1"}; !slices.Equal(bx, want) {
		t.Errorf("BAM output barcodes = %v, want %v", bx, want)
	}
}

func TestReport(t *testing.T) {
	var r report
	for _, splits := range []int{0, 0, 3, 1, 0, 3, 10} {
		r.add(splits)
	}
	var sb strings.Builder
	if err := r.write(&sb); err != nil {
		t.Fatal(err)
	}
	// 7 barcodes; 0+0+3+1+0+3+10 = 17 splits, so 24 molecules afterwards; rows in order of splits
	want := "Original_Molecules\t7\nNew_Molecules\t24\nDeconvolutions\tBarcodes\n0\t3\n1\t1\n3\t2\n10\t1\n"
	if sb.String() != want {
		t.Errorf("got:\n%q\nwant:\n%q", sb.String(), want)
	}
}

// A molecule is closed only when the scan is MORE than the distance past its last read; at exactly the
// distance it is still open.
func TestDeconvolve_ClosedOnlyBeyondTheDistance(t *testing.T) {
	path := writeSam(t,
		read("a", "c1", 100, aa),                    // AA, ends at 109
		read("b", "c1", 210, "BX:Z:BB\tVX:i:1"),     // starts at 209, exactly 100 past: AA is still open
		samLine("sup1", 2048, "c1", 211, "10M", aa), // so this still belongs to AA
		read("c", "c1", 212, "BX:Z:BB\tVX:i:1"),     // starts at 211, 102 past: AA is now closed
		samLine("sup2", 2048, "c1", 213, "10M", aa), // so this belongs to the next molecule, AA-1
	)
	got, _, err := run(t, path, 100)
	if err != nil {
		t.Fatal(err)
	}
	if want := []string{"a=AA", "b=BB", "sup1=AA", "c=BB", "sup2=AA-1"}; !slices.Equal(got, want) {
		t.Errorf("barcodes = %v, want %v", got, want)
	}
}
