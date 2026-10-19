package stats

import (
	"container/heap"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ── golden outputs from the Python original ───────────────────────────────────
//
// testdata/small.d*.tsv were produced by running harpy's bx_stats_sam.py on testdata/small.sam
// (see testdata/regen.sh). Output order is not part of the contract, so both sides are compared sorted.

var goldenDistances = []int{0, 1, 300, 2000, 40000}

func goldenCheck(t *testing.T, path string, d int) {
	t.Helper()
	got, err := captureStdout(t, func() error { return StatsXam(path, d, 2) })
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
	for _, d := range goldenDistances {
		goldenCheck(t, "testdata/small.sam", d)
	}
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
	for _, d := range []int{0, 300} {
		goldenCheck(t, path, d)
	}
}

// ── molecule assembly (expected lines from the Python version) ────────────────

func TestStatsXam_SplitSuffixAndSupplementary(t *testing.T) {
	const bc = "BX:Z:AA\tVX:i:1"
	path := writeSam(t,
		samLine("a", "c1", 100, bc), // 99..109
		// a supplementary alignment credits its 7 bp and insert to the open molecule, nothing else
		samLineCigar("s", 2048, "c1", 105, "7M", 7, bc),
		samLine("b", "c1", 120, bc),  // gap 10, joins the molecule
		samLine("c", "c1", 5000, bc), // gap far beyond 100, starts AA-1
	)
	got := statsOf(t, path, 100)
	want := []string{
		"c1\tAA\t2\t99\t129\t30\t27\t27\t0.9\t0.9",
		"c1\tAA-1\t1\t4999\t5009\t10\t10\t10\t1.0\t1.0",
	}
	if strings.Join(sortedLines(strings.Join(got, "\n")), "\n") != strings.Join(sortedLines(strings.Join(want, "\n")), "\n") {
		t.Errorf("got:\n%s\nwant:\n%s", strings.Join(got, "\n"), strings.Join(want, "\n"))
	}
}

// After a long read and a short one inside it, a read 140 bp past the short one starts a new molecule
// even though it is only 100 bp past the long one: the gap is measured from the previous read.
func TestStatsXam_GapMeasuredFromPreviousRead(t *testing.T) {
	const bc = "BX:Z:AA\tVX:i:1"
	path := writeSam(t,
		samLineCigar("a", 0, "c1", 100, "100M", 100, bc),
		samLineCigar("b", 0, "c1", 150, "10M", 10, bc),
		samLineCigar("c", 0, "c1", 300, "10M", 10, bc),
	)
	got := statsOf(t, path, 100)
	want := []string{
		"c1\tAA\t2\t99\t199\t100\t110\t110\t1.0\t1.0",
		"c1\tAA-1\t1\t299\t309\t10\t10\t10\t1.0\t1.0",
	}
	if strings.Join(sortedLines(strings.Join(got, "\n")), "\n") != strings.Join(sortedLines(strings.Join(want, "\n")), "\n") {
		t.Errorf("got:\n%s\nwant:\n%s", strings.Join(got, "\n"), strings.Join(want, "\n"))
	}
}

func TestStatsXam_NotSorted(t *testing.T) {
	const bc = "BX:Z:AA\tVX:i:1"
	path := writeSam(t, samLine("a", "c1", 100, bc), samLine("b", "c1", 50, bc))
	out, err := captureStdout(t, func() error { return StatsXam(path, 0, 1) })
	if !errors.Is(err, ErrNotSorted) {
		t.Errorf("err = %v, want ErrNotSorted", err)
	}
	if out != header {
		t.Errorf("output before the error should be just the header, got %q", out)
	}

	// a new contig may start at a lower position
	path = writeSam(t, samLine("a", "c1", 500, bc), samLine("b", "c2", 10, bc))
	if lines := statsOf(t, path, 0); len(lines) != 2 {
		t.Errorf("two contigs should give two molecules, got %q", lines)
	}
}

func TestStatsXam_EmptyFile(t *testing.T) {
	out, err := captureStdout(t, func() error { return StatsXam(writeSam(t), 0, 1) })
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
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path, 0), "\n")), "\n")
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
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path, 0), "\n")), "\n")
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
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path, 0), "\n")), "\n")
	want := strings.Join(sortedLines(strings.Join([]string{
		"c1\tACGTACGTACGTACGTAC\t2\t99\t209\t110\t20\t20\t0.18182\t0.18182",
		"c1\t1_2_3\t2\t299\t409\t110\t20\t20\t0.18182\t0.18182",
		"c1\tinvalid\t3\t0\t0\t0\t30\t0\t0.0\t0.0",
	}, "\n")), "\n")
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// The deadline list must hand back the entry furthest to the left first, and break ties by barcode name.
func TestEvictHeapOrder(t *testing.T) {
	var h evictHeap
	for _, e := range []evictEntry{{50, "b"}, {10, "z"}, {30, "a"}, {10, "a"}, {50, "a"}, {20, "m"}} {
		heap.Push(&h, e)
	}
	var got []string
	for h.Len() > 0 {
		e := heap.Pop(&h).(evictEntry)
		got = append(got, fmt.Sprintf("%d%s", e.end, e.bx))
	}
	if want := "10a 10z 20m 30a 50a 50b"; strings.Join(got, " ") != want {
		t.Errorf("pop order = %v, want %s", got, want)
	}
}

// A barcode that reappears after being closed early keeps counting its molecules (BARCODE-1, -2, ...).
func TestStatsXam_NumberingContinuesAfterEarlyClose(t *testing.T) {
	const bc = "BX:Z:AA\tVX:i:1"
	path := writeSam(t,
		samLine("a", "c1", 100, bc),   // AA: 99..109
		samLine("x", "c1", 5000, ""),  // moves the scan far past AA, so AA is closed early
		samLine("b", "c1", 5100, bc),  // AA again: molecule AA-1
		samLine("y", "c1", 10000, ""), // closes AA-1 early
		samLine("c", "c1", 10100, bc), // AA again: molecule AA-2
	)
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path, 100), "\n")), "\n")
	want := strings.Join(sortedLines(strings.Join([]string{
		"c1\tAA\t1\t99\t109\t10\t10\t10\t1.0\t1.0",
		"c1\tAA-1\t1\t5099\t5109\t10\t10\t10\t1.0\t1.0",
		"c1\tAA-2\t1\t10099\t10109\t10\t10\t10\t1.0\t1.0",
		"c1\tinvalid\t2\t0\t0\t0\t20\t0\t0.0\t0.0",
	}, "\n")), "\n")
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// A barcode that gets a newer read leaves its older line on the deadline list. When that old line becomes
// overdue it must be recognised as stale and discarded, not treated as the barcode being finished: here AA's
// first line (read ending at 109) is overdue when BB arrives at 249, but AA's newer read ended at 159 and the
// later AA read at 254 is still within 100 bp of it, so all three AA reads are one molecule.
func TestStatsXam_StaleDeadlineEntryIsIgnored(t *testing.T) {
	const aa, bb = "BX:Z:AA\tVX:i:1", "BX:Z:BB\tVX:i:1"
	path := writeSam(t,
		samLine("a1", "c1", 100, aa),
		samLine("a2", "c1", 150, aa),
		samLine("b", "c1", 250, bb),
		samLine("a3", "c1", 255, aa),
	)
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path, 100), "\n")), "\n")
	want := strings.Join(sortedLines(strings.Join([]string{
		"c1\tAA\t3\t99\t264\t165\t30\t30\t0.18182\t0.18182",
		"c1\tBB\t1\t249\t259\t10\t10\t10\t1.0\t1.0",
	}, "\n")), "\n")
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

// A gap of exactly the threshold does not split: only a gap bigger than it does. So when BB arrives at the
// very position where AA's gap equals the threshold, AA must not be closed yet.
func TestStatsXam_GapEqualToThresholdDoesNotSplit(t *testing.T) {
	const aa, bb = "BX:Z:AA\tVX:i:1", "BX:Z:BB\tVX:i:1"
	path := writeSam(t,
		samLine("a", "c1", 100, aa),
		samLine("b", "c1", 210, bb),
		samLine("a2", "c1", 210, aa),
	)
	got := strings.Join(sortedLines(strings.Join(statsOf(t, path, 100), "\n")), "\n")
	want := strings.Join(sortedLines(strings.Join([]string{
		"c1\tAA\t2\t99\t219\t120\t20\t20\t0.16667\t0.16667",
		"c1\tBB\t1\t209\t219\t10\t10\t10\t1.0\t1.0",
	}, "\n")), "\n")
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}
