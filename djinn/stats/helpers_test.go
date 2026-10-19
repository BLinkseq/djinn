package stats

import (
	"bytes"
	"testing"

	"github.com/biogo/hts/sam"
)

// The expected values below come from running pysam / Python, which the molecule table has to match.

// alignedLength / inferredQueryLength: cigar, SEQ length (0 = absent), then pysam's query_alignment_length and infer_query_length
var lengthCases = []struct {
	cigar                     string
	seqLen, aligned, inferred int
}{
	{"100M", 100, 100, 100},
	{"100M", 0, 100, 100},
	{"5S95M", 100, 95, 100},
	{"5S95M", 0, 95, 100},
	{"90M10S", 100, 90, 100},
	{"90M10S", 0, 90, 100},
	{"5S85M10S", 100, 85, 100},
	{"5S85M10S", 0, 85, 100},
	{"30M2I68M", 100, 100, 100},
	{"30M2I68M", 0, 100, 100},
	{"30M5D70M", 100, 100, 100},
	{"30M5D70M", 0, 100, 100},
	{"40M300N60M", 100, 100, 100},
	{"40M300N60M", 0, 100, 100},
	{"5H10S85M", 95, 85, 95},
	{"5H10S85M", 0, 85, 95},
	{"85M10S5H", 95, 85, 95},
	{"85M10S5H", 0, 85, 95},
	{"5H95M", 95, 95, 95},
	{"5H95M", 0, 95, 95},
	{"50=2X48=", 100, 100, 100},
	{"50=2X48=", 0, 100, 100},
	{"18M1P41M", 59, 59, 59},
	{"18M1P41M", 0, 59, 59},
	{"6S60M4D4M5S", 75, 64, 75},
	{"6S60M4D4M5S", 0, 64, 75},
	{"52M25S", 77, 52, 77},
	{"52M25S", 0, 52, 77},
	{"5S94M1S", 100, 94, 100},
	{"5S94M1S", 0, 94, 100},
	{"2M14S19H", 16, 2, 16},
	{"2M14S19H", 0, 2, 16},
	{"4I20M", 24, 24, 24},
	{"4I20M", 0, 24, 24},
	{"20M4I", 24, 24, 24},
	{"20M4I", 0, 24, 24},
	{"8H3S40M2S4H", 45, 40, 45},
	{"8H3S40M2S4H", 0, 40, 45},
	{"60M", 60, 60, 60},
	{"60M", 0, 60, 60},
	{"3S2I50M", 55, 52, 55},
	{"3S2I50M", 0, 52, 55},
}

// coverage: the Python expression max(0.0, round(min(num/denom, 1.0), 5)) printed with repr
var coverageCases = []struct {
	num, denom int
	want       string
}{
	{0, 5, "0.0"},
	{1, 3, "0.33333"},
	{2, 3, "0.66667"},
	{1, 7, "0.14286"},
	{5, 100000, "5e-05"},
	{1, 200000, "1e-05"},
	{3, 200000, "2e-05"},
	{1, 1, "1.0"},
	{7, 5, "1.0"},
	{123455, 1000000, "0.12345"},
	{123445, 1000000, "0.12344"},
	{1, 1000000, "0.0"},
	{4, 100000, "4e-05"},
	{5, 99999, "5e-05"},
	{99999, 100000, "0.99999"},
	{99999, 100001, "0.99998"},
	{1, 16, "0.0625"},
	{1, 32, "0.03125"},
	{3, 64, "0.04688"},
	{12, 13, "0.92308"},
	{1, 8, "0.125"},
	{7, 100000, "7e-05"},
	{2652, 632708, "0.00419"},
	{3234, 2730218, "0.00118"},
	{395, 303820, "0.0013"},
	{4389, 394811, "0.01112"},
	{2995, 2444391, "0.00123"},
	{475, 2128340, "0.00022"},
	{1758, 157269, "0.01118"},
	{704, 1818842, "0.00039"},
	{3425, 292995, "0.01169"},
	{1971, 380478, "0.00518"},
	{4514, 1780563, "0.00254"},
	{484, 2371685, "0.0002"},
	{1014, 936333, "0.00108"},
	{4775, 259469, "0.0184"},
	{4727, 2455939, "0.00192"},
	{3249, 207993, "0.01562"},
	{1811, 195382, "0.00927"},
	{4560, 558576, "0.00816"},
	{2372, 1757997, "0.00135"},
	{1181, 2267801, "0.00052"},
	{964, 2394586, "0.0004"},
	{2527, 2349890, "0.00108"},
	{1480, 432247, "0.00342"},
	{4764, 2395805, "0.00199"},
	{1539, 1561949, "0.00099"},
	{798, 2297407, "0.00035"},
	{514, 2367133, "0.00022"},
	{488, 2596315, "0.00019"},
	{1687, 2082114, "0.00081"},
	{4355, 1793453, "0.00243"},
	{2573, 1952876, "0.00132"},
	{4796, 1900794, "0.00252"},
	{2962, 1257314, "0.00236"},
	{2035, 753997, "0.0027"},
	{1999, 343325, "0.00582"},
	{4705, 1259337, "0.00374"},
	{4302, 2076669, "0.00207"},
	{2813, 1882548, "0.00149"},
	{2358, 2554160, "0.00092"},
	{599, 495204, "0.00121"},
}

func record(t *testing.T, cigar string, seqLen int) *sam.Record {
	t.Helper()
	c, err := sam.ParseCigar([]byte(cigar))
	if err != nil {
		t.Fatal(err)
	}
	r := &sam.Record{Cigar: c}
	if seqLen > 0 {
		r.Seq = sam.NewSeq(bytes.Repeat([]byte("A"), seqLen))
	}
	return r
}

func TestLengthsMatchPysam(t *testing.T) {
	for _, c := range lengthCases {
		r := record(t, c.cigar, c.seqLen)
		if got := alignedLength(r); got != c.aligned {
			t.Errorf("alignedLength(%s, seq %d) = %d, pysam says %d", c.cigar, c.seqLen, got, c.aligned)
		}
		if got := inferredQueryLength(r); got != c.inferred {
			t.Errorf("inferredQueryLength(%s) = %d, pysam says %d", c.cigar, got, c.inferred)
		}
	}
}

func TestInsertSize(t *testing.T) {
	r := record(t, "5S90M5S", 100)
	for _, c := range []struct {
		name  string
		flags sam.Flags
		tlen  int
		want  int
	}{
		{"unpaired, template longer", 0, -300, 300},
		{"unpaired, template shorter than the read", 0, 20, 100},
		{"unpaired, no template", 0, 0, 100},
		{"paired takes the template length", sam.Paired, 250, 250},
		{"paired negative template is clamped", sam.Paired, -250, 0},
		{"paired supplementary uses the aligned length", sam.Paired | sam.Supplementary, 250, 90},
		{"unpaired supplementary still uses the template", sam.Supplementary, -300, 300},
	} {
		r.Flags, r.TempLen = c.flags, c.tlen
		if got := insertSize(r); got != c.want {
			t.Errorf("%s: got %d, want %d", c.name, got, c.want)
		}
	}
}

func TestPyFloat(t *testing.T) {
	for f, want := range map[float64]string{
		0: "0.0", 1: "1.0", 0.5: "0.5", 0.99999: "0.99999", 0.0001: "0.0001", 0.00012: "0.00012",
		0.00001: "1e-05", 0.00003: "3e-05", 0.000012: "1.2e-05", 0.1: "0.1", 0.12346: "0.12346",
	} {
		if got := pyFloat(f); got != want {
			t.Errorf("pyFloat(%v) = %q, want %q", f, got, want)
		}
	}
}

func TestCoverageMatchesPython(t *testing.T) {
	for _, c := range coverageCases {
		if got := pyFloat(coverage(c.num, c.denom)); got != c.want {
			t.Errorf("coverage(%d, %d) = %s, Python says %s", c.num, c.denom, got, c.want)
		}
	}
}
