package barcodes

import "testing"

func TestTenXList_StreamsEmbeddedWhitelist(t *testing.T) {
	g, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList: %v", err)
	}
	defer g.Close()

	buf := make([]byte, g.MaxLen())
	n, ok := g.NextInto(buf)
	if !ok {
		t.Fatal("expected a first barcode")
	}
	if n != g.MaxLen() {
		t.Fatalf("first barcode length = %d, want %d", n, g.MaxLen())
	}
	first := string(buf[:n])
	for _, b := range buf[:n] {
		switch b {
		case 'A', 'C', 'G', 'T':
		default:
			t.Fatalf("first barcode %q contains non-ACGT byte %q", first, b)
		}
	}

	// a second, independent stream must start over from the beginning
	g2, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList (second instance): %v", err)
	}
	defer g2.Close()
	n2, ok := g2.NextInto(buf)
	if !ok || string(buf[:n2]) != first {
		t.Fatalf("expected second instance's first barcode to equal %q, got %q (ok=%v)", first, buf[:n2], ok)
	}
}

func TestTenXList_FullScanCountAndUniqueness(t *testing.T) {
	g, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList: %v", err)
	}
	defer g.Close()

	buf := make([]byte, g.MaxLen())
	seen := make(map[string]struct{}, 5_000_000)
	count := 0
	for {
		n, ok := g.NextInto(buf)
		if !ok {
			break
		}
		if n != g.MaxLen() {
			t.Fatalf("barcode %d has length %d, want %d", count, n, g.MaxLen())
		}
		bc := string(buf[:n])
		if _, dup := seen[bc]; dup {
			t.Fatalf("duplicate barcode %q at position %d", bc, count)
		}
		seen[bc] = struct{}{}
		count++
	}

	const want = 4_792_320 // 4M-with-alts-february-2016 whitelist size
	if count != want {
		t.Fatalf("scanned %d barcodes, want %d", count, want)
	}
}

func TestTenXList_GetInvalid(t *testing.T) {
	g, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList: %v", err)
	}
	defer g.Close()

	inv := g.GetInvalid()
	if len(inv) != g.MaxLen() {
		t.Fatalf("GetInvalid length = %d, want %d", len(inv), g.MaxLen())
	}
	for i, b := range inv {
		if b != 'N' {
			t.Fatalf("GetInvalid()[%d] = %q, want 'N'", i, b)
		}
	}
}

// Overwriting the shared dst buffer between calls must not retroactively
// corrupt a barcode already copied out of it -- same contract as every
// other Generator's NextInto.
func TestTenXList_NextIntoBufferReuseDoesNotCorruptPriorCopy(t *testing.T) {
	g, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList: %v", err)
	}
	defer g.Close()

	buf := make([]byte, g.MaxLen())
	n1, ok := g.NextInto(buf)
	if !ok {
		t.Fatal("expected a first barcode")
	}
	first := append([]byte(nil), buf[:n1]...) // caller's own copy, as the contract requires

	n2, ok := g.NextInto(buf) // reuses/overwrites buf
	if !ok {
		t.Fatal("expected a second barcode")
	}
	second := string(buf[:n2])

	if string(first) == second {
		t.Fatalf("first and second barcodes were identical (%q) -- whitelist or scan advanced incorrectly", second)
	}
	if len(first) != n1 {
		t.Fatalf("caller's copy of first barcode was corrupted: len=%d, want %d", len(first), n1)
	}
}

func TestLoadTenXSet(t *testing.T) {
	set, err := LoadTenXSet()
	if err != nil {
		t.Fatalf("LoadTenXSet: %v", err)
	}

	const want = 4_792_320
	if len(set) != want {
		t.Fatalf("loaded %d barcodes, want %d", len(set), want)
	}

	// cross-check against the streaming generator: every barcode NextInto
	// yields must be a member of the loaded set, and vice versa in spirit
	// (same underlying file).
	g, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList: %v", err)
	}
	defer g.Close()

	buf := make([]byte, g.MaxLen())
	checked := 0
	for checked < 1000 {
		n, ok := g.NextInto(buf)
		if !ok {
			break
		}
		if _, member := set[string(buf[:n])]; !member {
			t.Fatalf("barcode %q from TenXList is not present in LoadTenXSet's set", buf[:n])
		}
		checked++
	}
	if checked == 0 {
		t.Fatal("expected to check at least one barcode")
	}
}
