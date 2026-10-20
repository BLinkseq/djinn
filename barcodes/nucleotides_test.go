package barcodes

import "testing"

func TestNewGenericBounds(t *testing.T) {
	if _, err := NewGeneric(0); err == nil {
		t.Errorf("NewGeneric(0) expected error, got nil")
	}
	if _, err := NewGeneric(-1); err == nil {
		t.Errorf("NewGeneric(-1) expected error, got nil")
	}
	if _, err := NewGeneric(nucMaxLen + 1); err == nil {
		t.Errorf("NewGeneric(nucMaxLen+1) expected error, got nil")
	}
	if _, err := NewGeneric(nucMaxLen); err != nil {
		t.Errorf("NewGeneric(nucMaxLen) unexpected error: %v", err)
	}
	if _, err := NewGeneric(1); err != nil {
		t.Errorf("NewGeneric(1) unexpected error: %v", err)
	}
}

func TestNucleotidesNextLenOne(t *testing.T) {
	n, err := NewGeneric(1)
	if err != nil {
		t.Fatalf("NewGeneric(1): %v", err)
	}
	defer n.Close()

	want := []string{"A", "T", "C", "G"}
	for i, w := range want {
		bc, ok := n.Next()
		if !ok {
			t.Fatalf("expected barcode %d to be produced", i)
		}
		if got := bc.String(); got != w {
			t.Errorf("barcode %d = %q, want %q", i, got, w)
		}
	}
	// exhausted after 4 combinations
	if _, ok := n.Next(); ok {
		t.Errorf("expected sequence to be exhausted after 4 barcodes of length 1")
	}
}

func TestNucleotidesNextLenTwo(t *testing.T) {
	n, err := NewGeneric(2)
	if err != nil {
		t.Fatalf("NewGeneric(2): %v", err)
	}
	defer n.Close()

	// alphabet order is A,T,C,G so length-2 sequence starts AA, AT, AC, AG, TA...
	want := []string{"AA", "AT", "AC", "AG", "TA"}
	for i, w := range want {
		bc, ok := n.Next()
		if !ok {
			t.Fatalf("expected barcode %d to be produced", i)
		}
		if got := bc.String(); got != w {
			t.Errorf("barcode %d = %q, want %q", i, got, w)
		}
	}
}

func TestNucleotidesTotalCount(t *testing.T) {
	n, err := NewGeneric(2)
	if err != nil {
		t.Fatalf("NewGeneric(2): %v", err)
	}
	defer n.Close()

	count := 0
	seen := map[string]bool{}
	for {
		bc, ok := n.Next()
		if !ok {
			break
		}
		seen[bc.String()] = true
		count++
		if count > 100 {
			t.Fatalf("sequence did not terminate as expected")
		}
	}
	if count != 16 { // 4^2
		t.Errorf("got %d barcodes, want 16", count)
	}
	if len(seen) != 16 {
		t.Errorf("got %d distinct barcodes, want 16 (duplicates present)", len(seen))
	}
}

func TestNucleotidesGetInvalid(t *testing.T) {
	n, err := NewGeneric(4)
	if err != nil {
		t.Fatalf("NewGeneric(4): %v", err)
	}
	defer n.Close()
	if got, want := string(n.GetInvalid()), "NNNN"; got != want {
		t.Errorf("GetInvalid() = %q, want %q", got, want)
	}
}

func TestNewTellseqLength(t *testing.T) {
	n := NewTellseq()
	defer n.Close()
	if got, want := n.MaxLen(), 18; got != want {
		t.Errorf("MaxLen() = %d, want %d", got, want)
	}
	if got, want := len(n.GetInvalid()), 18; got != want {
		t.Errorf("len(GetInvalid()) = %d, want %d", got, want)
	}
	bc, ok := n.Next()
	if !ok {
		t.Fatalf("expected barcode to be produced")
	}
	if got, want := len(bc.Bytes()), 18; got != want {
		t.Errorf("first barcode length = %d, want %d", got, want)
	}
}

func TestNewTenXLength(t *testing.T) {
	n := NewTenX()
	defer n.Close()
	if got, want := n.MaxLen(), 16; got != want {
		t.Errorf("MaxLen() = %d, want %d", got, want)
	}
	if got, want := len(n.GetInvalid()), 16; got != want {
		t.Errorf("len(GetInvalid()) = %d, want %d", got, want)
	}
}

func TestNucleotidesNextInto(t *testing.T) {
	n, err := NewGeneric(3)
	if err != nil {
		t.Fatalf("NewGeneric(3): %v", err)
	}
	defer n.Close()
	dst := make([]byte, n.MaxLen())
	written, ok := n.NextInto(dst)
	if !ok {
		t.Fatalf("expected barcode to be produced")
	}
	if got, want := string(dst[:written]), "AAA"; got != want {
		t.Errorf("NextInto wrote %q, want %q", got, want)
	}
}
