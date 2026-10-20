package barcodes

import "testing"

func TestStlfrNext(t *testing.T) {
	s := NewStlfr()
	defer s.Close()

	bc, ok := s.Next()
	if !ok {
		t.Fatalf("expected first barcode to be produced")
	}
	if got, want := bc.String(), "1_1_1"; got != want {
		t.Errorf("first barcode = %q, want %q", got, want)
	}

	bc, ok = s.Next()
	if !ok {
		t.Fatalf("expected second barcode to be produced")
	}
	if got, want := bc.String(), "1_1_2"; got != want {
		t.Errorf("second barcode = %q, want %q", got, want)
	}
}

func TestStlfrBytesAndString(t *testing.T) {
	bc := Barcode{Data: [stlfrMaxLen]byte{'1', '_', '2', '_', '3'}, Len: 5}
	if got, want := bc.String(), "1_2_3"; got != want {
		t.Errorf("String() = %q, want %q", got, want)
	}
	if got, want := string(bc.Bytes()), "1_2_3"; got != want {
		t.Errorf("Bytes() = %q, want %q", got, want)
	}
}

func TestStlfrGetInvalid(t *testing.T) {
	s := NewStlfr()
	defer s.Close()
	if got, want := string(s.GetInvalid()), "0_0_0"; got != want {
		t.Errorf("GetInvalid() = %q, want %q", got, want)
	}
}

func TestStlfrInvalidInto(t *testing.T) {
	s := NewStlfr()
	defer s.Close()
	dst := make([]byte, s.MaxLen())
	n := s.InvalidInto(dst)
	if got, want := string(dst[:n]), "0_0_0"; got != want {
		t.Errorf("InvalidInto wrote %q, want %q", got, want)
	}
}

func TestStlfrNextInto(t *testing.T) {
	s := NewStlfr()
	defer s.Close()
	dst := make([]byte, s.MaxLen())
	n, ok := s.NextInto(dst)
	if !ok {
		t.Fatalf("expected barcode to be produced")
	}
	if got, want := string(dst[:n]), "1_1_1"; got != want {
		t.Errorf("NextInto wrote %q, want %q", got, want)
	}
}

func TestStlfrMaxLen(t *testing.T) {
	s := NewStlfr()
	defer s.Close()
	if got, want := s.MaxLen(), 14; got != want {
		t.Errorf("MaxLen() = %d, want %d", got, want)
	}
}

// TestStlfrExhaustion walks the sequence until the last value of the
// innermost dimension rolls over into the next value of the middle
// dimension, verifying the iteration order and bounds (1..1537 inclusive)
// without iterating the full 1537^3 space.
func TestStlfrExhaustion(t *testing.T) {
	s := NewStlfr()
	defer s.Close()

	// advance to the end of the innermost "z" loop for x=1,y=1
	var last Barcode
	for i := 0; i < 1537; i++ {
		bc, ok := s.Next()
		if !ok {
			t.Fatalf("sequence ended prematurely at i=%d", i)
		}
		last = bc
	}
	if got, want := last.String(), "1_1_1537"; got != want {
		t.Errorf("1537th barcode = %q, want %q", got, want)
	}

	// next call should roll over to y=2
	bc, ok := s.Next()
	if !ok {
		t.Fatalf("expected barcode after rollover")
	}
	if got, want := bc.String(), "1_2_1"; got != want {
		t.Errorf("barcode after rollover = %q, want %q", got, want)
	}
}

func TestStlfrCloseStopsIteration(t *testing.T) {
	s := NewStlfr()
	s.Close()
	// Calling Next after Close should not panic; iter.Pull's stop function
	// makes subsequent next() calls return the zero value and false.
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("Next() after Close panicked: %v", r)
		}
	}()
	_, _ = s.Next()
}
