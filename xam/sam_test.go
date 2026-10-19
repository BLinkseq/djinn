package xam

import (
	"testing"

	"github.com/biogo/hts/sam"
)

func TestIsValid(t *testing.T) {
	cases := []struct {
		name    string
		barcode string
		want    bool
	}{
		{"clean haplotag", "A01C02B03D04", true},
		{"contains N", "A01C02B03DN4", false},
		{"has A00 segment", "A00C02B03D04", false},
		{"has B00 segment", "A01B00B03D04", false},
		{"has C00 segment", "A01C00B03D04", false},
		{"has D00 segment", "A01C02B03D00", false},
		{"stlfr prefix zero", "0_1_2", false},
		{"stlfr suffix zero", "1_2_0", false},
		{"stlfr middle zero", "1_0_2", false},
		{"stlfr valid", "1_2_3", true},
		{"empty string valid", "", true},
		{"tellseq clean", "ATCGATCG", true},
		{"tellseq with N", "ATCGNATCG", false},
		{"clean haplotag with hyphen suffix", "A01C02B03D04-1", true},
		{"haplotag A00 despite hyphen suffix", "A00C02B03D04-1", false},
		{"suffix content never affects validity", "1_2_3-NN", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := IsValid(c.barcode); got != c.want {
				t.Errorf("IsValid(%q) = %v, want %v", c.barcode, got, c.want)
			}
		})
	}
}

func TestBoolToInt(t *testing.T) {
	if BoolToInt(true) != 1 {
		t.Errorf("BoolToInt(true) = %d, want 1", BoolToInt(true))
	}
	if BoolToInt(false) != 0 {
		t.Errorf("BoolToInt(false) = %d, want 0", BoolToInt(false))
	}
}

func TestBoolToSInt(t *testing.T) {
	if got := BoolToSInt(true); got != "1" {
		t.Errorf("BoolToSInt(true) = %q, want %q", got, "1")
	}
	if got := BoolToSInt(false); got != "0" {
		t.Errorf("BoolToSInt(false) = %q, want %q", got, "0")
	}
}

func newRecord(t *testing.T, name string) *sam.Record {
	t.Helper()
	rec, err := sam.NewRecord(name, nil, nil, -1, -1, 0, 0, nil, []byte("ACGT"), nil, nil)
	if err != nil {
		t.Fatalf("sam.NewRecord: %v", err)
	}
	return rec
}

func TestSetBxStringAndGetStringTag(t *testing.T) {
	rec := newRecord(t, "read1")
	val := "ATCGATCG"
	SetBxString(rec, &val)

	got, ok := GetStringTag(rec, "BX")
	if !ok {
		t.Fatal("expected BX tag to be present")
	}
	if got != val {
		t.Errorf("GetStringTag = %q, want %q", got, val)
	}

	// setting again should replace, not duplicate
	val2 := "GGGGCCCC"
	SetBxString(rec, &val2)
	if len(rec.AuxFields) != 1 {
		t.Fatalf("expected 1 aux field after replace, got %d", len(rec.AuxFields))
	}
	got2, ok := GetStringTag(rec, "BX")
	if !ok || got2 != val2 {
		t.Errorf("GetStringTag after replace = %q, %v, want %q, true", got2, ok, val2)
	}
}

func TestSetBxByteAndGetStringTag(t *testing.T) {
	rec := newRecord(t, "read1")
	val := []byte("TTTTAAAA")
	SetBxByte(rec, &val)

	got, ok := GetStringTag(rec, "BX")
	if !ok {
		t.Fatal("expected BX tag to be present")
	}
	if got != string(val) {
		t.Errorf("GetStringTag = %q, want %q", got, val)
	}

	// replacing via SetBxByte should not duplicate the field
	val2 := []byte("CCCCGGGG")
	SetBxByte(rec, &val2)
	if len(rec.AuxFields) != 1 {
		t.Fatalf("expected 1 aux field after replace, got %d", len(rec.AuxFields))
	}
}

func TestGetStringTagMissing(t *testing.T) {
	rec := newRecord(t, "read1")
	got, ok := GetStringTag(rec, "BX")
	if ok || got != "" {
		t.Errorf("GetStringTag on empty record = %q, %v, want \"\", false", got, ok)
	}
}

func TestGetIntTagWidths(t *testing.T) {
	// GetIntTag must work regardless of which concrete width the aux
	// library chose to store the value as ('c'/'C'/'s'/'S'/'i'/'I').
	cases := []struct {
		name string
		aux  sam.Aux
		want int
	}{
		{"int8 'c'", sam.Aux{'V', 'X', 'c', 1}, 1},
		{"uint8 'C'", sam.Aux{'V', 'X', 'C', 1}, 1},
		{"int32 'i' zero", func() sam.Aux { a := sam.Aux{'V', 'X', 'i', 0, 0, 0, 0}; return a }(), 0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			rec := newRecord(t, "read1")
			rec.AuxFields = append(rec.AuxFields, c.aux)
			got, ok := GetIntTag(rec, "VX")
			if !ok {
				t.Fatalf("GetIntTag did not find tag for width %s", c.name)
			}
			if got != c.want {
				t.Errorf("GetIntTag = %d, want %d", got, c.want)
			}
		})
	}
}

func TestGetIntTagMissing(t *testing.T) {
	rec := newRecord(t, "read1")
	got, ok := GetIntTag(rec, "VX")
	if ok || got != 0 {
		t.Errorf("GetIntTag on empty record = %d, %v, want 0, false", got, ok)
	}
}

func TestSetVXAndGetVX(t *testing.T) {
	rec := newRecord(t, "read1")
	SetVX(rec, true)
	val, ok := GetVX(rec)
	if !ok {
		t.Fatal("expected VX tag present after SetVX")
	}
	if !val {
		t.Errorf("GetVX = false, want true")
	}

	SetVX(rec, false)
	if len(rec.AuxFields) != 1 {
		t.Fatalf("expected 1 aux field after re-set, got %d", len(rec.AuxFields))
	}
	val, ok = GetVX(rec)
	if !ok {
		t.Fatal("expected VX tag present after second SetVX")
	}
	if val {
		t.Errorf("GetVX = true, want false")
	}
}

func TestSetVXOnNarrowExistingField(t *testing.T) {
	// Regression test: an aux field parsed from SAM text such as "VX:i:1"
	// can be encoded internally with a narrower type (e.g. 'C', 1 byte)
	// rather than the 4-byte 'i' SetVX used to assume. Previously this
	// caused SetVX to panic with a slice-bounds-out-of-range error.
	rec := newRecord(t, "read1")
	rec.AuxFields = append(rec.AuxFields, sam.Aux{'V', 'X', 'C', 1})

	SetVX(rec, false)

	val, ok := GetVX(rec)
	if !ok {
		t.Fatal("expected VX tag present after SetVX on narrow field")
	}
	if val {
		t.Errorf("GetVX = true, want false")
	}
	if len(rec.AuxFields) != 1 {
		t.Errorf("expected exactly 1 aux field, got %d", len(rec.AuxFields))
	}
}

func TestGetVXNotPresent(t *testing.T) {
	rec := newRecord(t, "read1")
	val, ok := GetVX(rec)
	if ok {
		t.Fatal("expected VX not present")
	}
	if val {
		t.Errorf("GetVX value = true, want false (zero value)")
	}
}

func TestFindBarcodeFromBXTag(t *testing.T) {
	rec := newRecord(t, "read1")
	bx := "A01C02B03D04"
	SetBxString(rec, &bx)

	gotBX, gotVX := FindBarcode(rec)
	if gotBX != bx {
		t.Errorf("FindBarcode bx = %q, want %q", gotBX, bx)
	}
	if !gotVX {
		t.Errorf("FindBarcode vx = false, want true (valid barcode, no VX tag => inferred)")
	}
}

func TestFindBarcodeInfersInvalidFromBX(t *testing.T) {
	rec := newRecord(t, "read1")
	bx := "A00C02B03D04" // invalid: has A00 segment
	SetBxString(rec, &bx)

	gotBX, gotVX := FindBarcode(rec)
	if gotBX != bx {
		t.Errorf("FindBarcode bx = %q, want %q", gotBX, bx)
	}
	if gotVX {
		t.Errorf("FindBarcode vx = true, want false")
	}
}

func TestFindBarcodeFromNameTellseq(t *testing.T) {
	rec := newRecord(t, "READ123:ATCGATCG")
	gotBX, gotVX := FindBarcode(rec)
	if gotBX != "ATCGATCG" {
		t.Errorf("FindBarcode bx = %q, want %q", gotBX, "ATCGATCG")
	}
	if !gotVX {
		t.Errorf("FindBarcode vx = false, want true")
	}
}

func TestFindBarcodeFromNameStlfr(t *testing.T) {
	rec := newRecord(t, "READ123#1_2_3")
	gotBX, gotVX := FindBarcode(rec)
	if gotBX != "1_2_3" {
		t.Errorf("FindBarcode bx = %q, want %q", gotBX, "1_2_3")
	}
	if !gotVX {
		t.Errorf("FindBarcode vx = false, want true")
	}
}

func TestFindBarcodeFromNameTellseqWithHyphenSuffix(t *testing.T) {
	rec := newRecord(t, "READ123:ATCGATCG-1")
	gotBX, gotVX := FindBarcode(rec)
	if gotBX != "ATCGATCG-1" {
		t.Errorf("FindBarcode bx = %q, want %q", gotBX, "ATCGATCG-1")
	}
	if !gotVX {
		t.Errorf("FindBarcode vx = false, want true")
	}
}

func TestFindBarcodeFromNameStlfrWithHyphenSuffix(t *testing.T) {
	rec := newRecord(t, "READ123#1_2_3-1")
	gotBX, gotVX := FindBarcode(rec)
	if gotBX != "1_2_3-1" {
		t.Errorf("FindBarcode bx = %q, want %q", gotBX, "1_2_3-1")
	}
	if !gotVX {
		t.Errorf("FindBarcode vx = false, want true")
	}
}

func TestFindBarcodeNoneFound(t *testing.T) {
	rec := newRecord(t, "READ123")
	gotBX, gotVX := FindBarcode(rec)
	if gotBX != "" {
		t.Errorf("FindBarcode bx = %q, want empty", gotBX)
	}
	if gotVX {
		t.Errorf("FindBarcode vx = true, want false")
	}
}

func TestFindBarcodeExplicitVXHonored(t *testing.T) {
	rec := newRecord(t, "read1")
	bx := "A00C02B03D04" // would be invalid if inferred
	SetBxString(rec, &bx)
	SetVX(rec, true) // explicitly mark valid, should override inference

	gotBX, gotVX := FindBarcode(rec)
	if gotBX != bx {
		t.Errorf("FindBarcode bx = %q, want %q", gotBX, bx)
	}
	if !gotVX {
		t.Errorf("FindBarcode vx = false, want true (explicit VX tag should be honored)")
	}
}

func TestPairedFlag(t *testing.T) {
	if got := PairedFlag(true); got != 77 {
		t.Errorf("PairedFlag(true) = %d, want 77", got)
	}
	if got := PairedFlag(false); got != 141 {
		t.Errorf("PairedFlag(false) = %d, want 141", got)
	}
}

func TestAlignedEnd(t *testing.T) {
	for cigar, want := range map[string]int{
		"100M": 110, "5S95M": 105, "30M5D70M": 115, "40M300N60M": 410, "50=2X48=": 110,
		"8H3S40M2S4H": 50, "10I": 11, // no reference consumed: htslib counts one base
	} {
		c, err := sam.ParseCigar([]byte(cigar))
		if err != nil {
			t.Fatal(err)
		}
		r := &sam.Record{Pos: 10, Cigar: c}
		if got := AlignedEnd(r); got != want {
			t.Errorf("AlignedEnd(%s at 10) = %d, want %d", cigar, got, want)
		}
	}
}
