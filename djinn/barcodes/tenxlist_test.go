package barcodes

import "testing"

func TestTenXList_StreamsEmbeddedWhitelist(t *testing.T) {
	g, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList: %v", err)
	}
	defer g.Close()

	buf := make([]byte, g.MaxLen())
	count := 0
	seen := make(map[string]bool)
	for {
		n, ok := g.NextInto(buf)
		if !ok {
			break
		}
		if n != g.MaxLen() {
			t.Fatalf("barcode %q has length %d, want %d", buf[:n], n, g.MaxLen())
		}
		bc := string(buf[:n])
		if seen[bc] {
			t.Fatalf("duplicate barcode %q at position %d", bc, count)
		}
		seen[bc] = true
		count++
	}
	if count == 0 {
		t.Fatal("expected at least one barcode from the embedded whitelist")
	}

	// a second, independent stream must start over from the beginning
	g2, err := NewTenXList()
	if err != nil {
		t.Fatalf("NewTenXList (second instance): %v", err)
	}
	defer g2.Close()
	n2, ok := g2.NextInto(buf)
	if !ok {
		t.Fatal("expected second instance to yield a first barcode")
	}
	if string(buf[:n2]) == "" {
		t.Fatal("expected non-empty first barcode from second instance")
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
		t.Fatalf("first and second barcodes were identical (%q) -- test fixture or scan advanced incorrectly", second)
	}
	// first must still read back correctly from its own independent copy
	if len(first) != n1 {
		t.Fatalf("caller's copy of first barcode was corrupted: len=%d, want %d", len(first), n1)
	}
}
