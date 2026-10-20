package barcodes

import "testing"

func TestWriteHaplotagPart(t *testing.T) {
	tests := []struct {
		n    int
		want string
	}{
		{1, "A01"},
		{9, "A09"},
		{10, "A10"},
		{96, "A96"},
	}
	for _, tt := range tests {
		var buf [12]byte
		writeHaplotagPart(&buf, 0, 'A', tt.n)
		if got := string(buf[0:3]); got != tt.want {
			t.Errorf("writeHaplotagPart(n=%d) = %q, want %q", tt.n, got, tt.want)
		}
	}
}

func TestHaplotaggingNext(t *testing.T) {
	h := NewHaplotagging()
	defer h.Close()

	bc, ok := h.Next()
	if !ok {
		t.Fatalf("expected first barcode to be produced")
	}
	if got, want := string(bc[:]), "A01C01B01D01"; got != want {
		t.Errorf("first barcode = %q, want %q", got, want)
	}

	bc, ok = h.Next()
	if !ok {
		t.Fatalf("expected second barcode to be produced")
	}
	if got, want := string(bc[:]), "A01C01B01D02"; got != want {
		t.Errorf("second barcode = %q, want %q", got, want)
	}
}

func TestHaplotaggingRollover(t *testing.T) {
	h := NewHaplotagging()
	defer h.Close()

	// advance through all 96 D values for A01C01B01
	var last [12]byte
	for i := 0; i < 96; i++ {
		bc, ok := h.Next()
		if !ok {
			t.Fatalf("sequence ended prematurely at i=%d", i)
		}
		last = bc
	}
	if got, want := string(last[:]), "A01C01B01D96"; got != want {
		t.Errorf("96th barcode = %q, want %q", got, want)
	}

	// next should roll D back to 01 and bump B to 02
	bc, ok := h.Next()
	if !ok {
		t.Fatalf("expected barcode after rollover")
	}
	if got, want := string(bc[:]), "A01C01B02D01"; got != want {
		t.Errorf("barcode after rollover = %q, want %q", got, want)
	}
}

func TestHaplotaggingGetInvalid(t *testing.T) {
	h := NewHaplotagging()
	defer h.Close()
	if got, want := string(h.GetInvalid()), "A00C00B00D00"; got != want {
		t.Errorf("GetInvalid() = %q, want %q", got, want)
	}
}

func TestHaplotaggingMaxLen(t *testing.T) {
	h := NewHaplotagging()
	defer h.Close()
	if got, want := h.MaxLen(), 12; got != want {
		t.Errorf("MaxLen() = %d, want %d", got, want)
	}
}

func TestHaplotaggingNextInto(t *testing.T) {
	h := NewHaplotagging()
	defer h.Close()
	dst := make([]byte, h.MaxLen())
	n, ok := h.NextInto(dst)
	if !ok {
		t.Fatalf("expected barcode to be produced")
	}
	if got, want := string(dst[:n]), "A01C01B01D01"; got != want {
		t.Errorf("NextInto wrote %q, want %q", got, want)
	}
}

func TestHaplotaggingNextIntoShortBuffer(t *testing.T) {
	h := NewHaplotagging()
	defer h.Close()
	dst := make([]byte, 5)
	n, ok := h.NextInto(dst)
	if !ok {
		t.Fatalf("expected barcode to be produced")
	}
	// copy() truncates to len(dst); verify it doesn't panic and reports
	// the truncated length actually written.
	if n != 5 {
		t.Errorf("NextInto into short buffer returned n=%d, want 5", n)
	}
}
