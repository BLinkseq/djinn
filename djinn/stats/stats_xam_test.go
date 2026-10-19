package stats

import (
	"errors"
	"fmt"
	"io"
	"os"
	"slices"
	"strings"
	"testing"

	"djinn/deconvolve"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ── golden outputs from the Python original ───────────────────────────────────
//
// testdata/small.d*.tsv were produced by running harpy's bx_stats_sam.py on testdata/small.sam
// (see testdata/regen.sh). Output order is not part of the contract, so both sides are compared sorted.

// goldenCheck runs StatsXam on path and compares it with the Python output for distance d
// (testdata/small.d<d>.tsv). Only the plain file (d = 0) is a direct match; the others are the Python
// version's output with its built-in splitting, which in djinn is the job of `sam deconvolve`
// (see TestStatsXam_AfterDeconvolve).
func goldenCheck(t *testing.T, path string, d int) {
	t.Helper()
	got, err := captureStdout(t, func() error { return StatsXam(path, 2) })
	if err != nil {
		t.Fatalf("d=%d: %v", d, err)
	}
	want, err := os.ReadFile(fmt.Sprintf("testdata/small.d%d.tsv", d))
	if err != nil {
		t.Fatal(err)
	}
	g, w := sortedLines(got), sortedLines(string(want))
	if strings.Join(g, "\n") != strings.Join(w, "\n") {
		t.Errorf("d=%d: output differs from the Python version (%d vs %d lines)", d, len(g), len(w))
		for i := 0; i < len(g) && i < len(w); i++ {
			if g[i] != w[i] {
				t.Logf("first difference:\n  go:     %q\n  python: %q", g[i], w[i])
				break
			}
		}
	}
}

func TestStatsXam_GoldenSAM(t *testing.T) {
	goldenCheck(t, "testdata/small.sam", 0)
}

func TestStatsXam_GoldenBAM(t *testing.T) {
	f, err := os.Open("testdata/small.sam")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	sr, err := sam.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	path := tempPath(t, "small.bam")
	out, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	bw, err := bam.NewWriter(out, sr.Header(), 1)
	if err != nil {
		t.Fatal(err)
	}
	for {
		rec, err := sr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		if err := bw.Write(rec); err != nil {
			t.Fatal(err)
		}
	}
	if err := errors.Join(bw.Close(), out.Close()); err != nil {
		t.Fatal(err)
	}
	goldenCheck(t, path, 0)
}

// ── molecule assembly ─────────────────────────────────────────────────────────

// Supplementary alignments credit their bases to the molecule but not its read count or extent (expected
// line from the Python version).
func TestStatsXam_Supplementary(t *testing.T) {
	const bc = "BX:Z:AA\tVX:i:1"
	path := writeSam(t,
		samLine("a", "c1", 100, bc), // 99..109
		samLineCigar("s", 2048, "c1", 105, "7M", 7, bc),
		samLine("b", "c1", 120, bc),
	)
	got := statsOf(t, path)
	want := []string{"c1\tAA\t2\t99\t129\t30\t27\t27\t0.9\t0.9"}
	if !slices.Equal(got, want) {
		t.Errorf("got:\n%s\nwant:\n%s", strings.Join(got, "\n"), strings.Join(want, "\n"))
	}
}

// Numbered barcodes (what `sam deconvolve` writes) are separate molecules.
func TestStatsXam_NumberedBarcodesAreSeparateMolecules(t *testing.T) {
	path := writeSam(t,
		samLine("a", "c1", 100, "BX:Z:AA\tVX:i:1"),
		samLine("b", "c1", 120, "BX:Z:AA\tVX:i:1"),
		samLine("c", "c1", 5000, "BX:Z:AA-1\tVX:i:1"),
		samLine("d", "c1", 5050, "BX:Z:AA-1"), // no VX: valid by the barcode itself, the "-1" does not matter
	)
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path), "\n")), "\n")
	want := strings.Join(sortedLines(strings.Join([]string{
		"c1\tAA\t2\t99\t129\t30\t20\t20\t0.66667\t0.66667",
		"c1\tAA-1\t2\t4999\t5059\t60\t20\t20\t0.33333\t0.33333",
	}, "\n")), "\n")
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// The split that used to be built into stats is now `sam deconvolve`. Deconvolving and then running stats on
// the result must give exactly what the Python version gave with its built-in -d.
func TestStatsXam_AfterDeconvolve(t *testing.T) {
	for _, d := range []int{1, 300, 2000, 40000} {
		// deconvolve the plain-text SAM, writing SAM so the next step can read it
		var sam string
		_, err := captureStderr(t, func() error {
			var err error
			sam, err = captureStdout(t, func() error { return deconvolve.DeconvolveXam("testdata/small.sam", d, true, 2) })
			return err
		})
		if err != nil {
			t.Fatalf("d=%d: deconvolve: %v", d, err)
		}
		path := tempPath(t, fmt.Sprintf("deconv.%d.sam", d))
		if err := os.WriteFile(path, []byte(sam), 0o644); err != nil {
			t.Fatal(err)
		}
		goldenCheck(t, path, d)
	}
}

func TestStatsXam_NotSorted(t *testing.T) {
	const bc = "BX:Z:AA\tVX:i:1"
	path := writeSam(t, samLine("a", "c1", 100, bc), samLine("b", "c1", 50, bc))
	out, err := captureStdout(t, func() error { return StatsXam(path, 1) })
	if !errors.Is(err, ErrNotSorted) {
		t.Errorf("err = %v, want ErrNotSorted", err)
	}
	if out != header {
		t.Errorf("output before the error should be just the header, got %q", out)
	}

	// a new contig may start at a lower position
	path = writeSam(t, samLine("a", "c1", 500, bc), samLine("b", "c2", 10, bc))
	if lines := statsOf(t, path); len(lines) != 2 {
		t.Errorf("two contigs should give two molecules, got %q", lines)
	}
}

func TestStatsXam_EmptyFile(t *testing.T) {
	out, err := captureStdout(t, func() error { return StatsXam(writeSam(t), 1) })
	if err != nil || out != header {
		t.Errorf("empty input: err = %v, output %q", err, out)
	}
}

func TestStatsXam_InvalidIsPooledPerContig(t *testing.T) {
	path := writeSam(t,
		samLine("a", "c1", 100, "BX:Z:AA\tVX:i:0"),
		samLine("b", "c1", 200, ""),
		samLine("c", "c2", 100, "BX:Z:AA\tVX:i:0"),
	)
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path), "\n")), "\n")
	want := "c1\tinvalid\t2\t0\t0\t0\t20\t0\t0.0\t0.0\nc2\tinvalid\t1\t0\t0\t0\t10\t0\t0.0\t0.0"
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// ── barcode identification through xam.FindBarcode ────────────────────────────

// With no VX tag a barcode's validity is decided by the barcode itself; a VX tag always wins.
func TestStatsXam_VXInferredFromBarcode(t *testing.T) {
	path := writeSam(t,
		samLine("r1", "c1", 100, "BX:Z:A01C02B03D04"), // no VX, valid barcode
		samLine("r2", "c1", 200, "BX:Z:A01C02B03D04"),
		samLine("r3", "c1", 300, "BX:Z:A00C02B03D04"),         // no VX, A00 -> invalid
		samLine("r4", "c1", 400, "BX:Z:ACGTNCGT"),             // no VX, contains N -> invalid
		samLine("r5", "c1", 500, "BX:Z:A09C02B03D04\tVX:i:1"), // explicit valid
		samLine("r6", "c1", 600, "BX:Z:A00C00B00D00\tVX:i:1"), // VX wins over an invalid-looking barcode
		samLine("r7", "c1", 700, "BX:Z:A01C01B01D01\tVX:i:0"), // VX wins over a valid-looking barcode
		samLine("r8", "c1", 800, ""),                          // no barcode at all
	)
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path), "\n")), "\n")
	want := strings.Join(sortedLines(strings.Join([]string{
		"c1\tA01C02B03D04\t2\t99\t209\t110\t20\t20\t0.18182\t0.18182",
		"c1\tA09C02B03D04\t1\t499\t509\t10\t10\t10\t1.0\t1.0",
		"c1\tA00C00B00D00\t1\t599\t609\t10\t10\t10\t1.0\t1.0",
		"c1\tinvalid\t4\t0\t0\t0\t40\t0\t0.0\t0.0", // r3, r4, r7, r8
	}, "\n")), "\n")
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// Without a BX tag the barcode is taken from the read name (TELLseq ":BASES" or stLFR "#n_n_n").
func TestStatsXam_BarcodeFromReadName(t *testing.T) {
	path := writeSam(t,
		samLine("r1:ACGTACGTACGTACGTAC", "c1", 100, ""),
		samLine("r2:ACGTACGTACGTACGTAC", "c1", 200, ""),
		samLine("r3#1_2_3", "c1", 300, ""),
		samLine("r4#1_2_3", "c1", 400, ""),
		samLine("r5#0_2_3", "c1", 500, ""),      // stLFR zero component -> invalid
		samLine("r6:ACGTNCGTAC", "c1", 600, ""), // contains N -> invalid
		samLine("plain", "c1", 700, ""),
	)
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path), "\n")), "\n")
	want := strings.Join(sortedLines(strings.Join([]string{
		"c1\tACGTACGTACGTACGTAC\t2\t99\t209\t110\t20\t20\t0.18182\t0.18182",
		"c1\t1_2_3\t2\t299\t409\t110\t20\t20\t0.18182\t0.18182",
		"c1\tinvalid\t3\t0\t0\t0\t30\t0\t0.0\t0.0",
	}, "\n")), "\n")
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}
