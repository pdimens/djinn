package fastq

import (
	"testing"

	"github.com/shenwei356/bio/seqio/fastx"
)

func newRecord(id, desc string) *fastx.Record {
	return &fastx.Record{
		ID:   []byte(id),
		Desc: []byte(desc),
	}
}

func TestHaplotagBX_ValidWithVX(t *testing.T) {
	rec := newRecord("read1", "VX:i:1\tBX:Z:A01C02B03D04")
	valid := HaplotagBX(rec)
	if !valid {
		t.Errorf("expected valid=true")
	}
	if got, want := string(rec.Desc), "VX:i:1\tBX:Z:A01C02B03D04"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestHaplotagBX_ValidNoVXInfersFromBarcode(t *testing.T) {
	rec := newRecord("read1", "BX:Z:A01C02B03D04")
	valid := HaplotagBX(rec)
	if !valid {
		t.Errorf("expected valid=true (no VX tag, barcode itself is valid)")
	}
	if got, want := string(rec.Desc), "VX:i:1\tBX:Z:A01C02B03D04"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestHaplotagBX_InvalidBarcodeNoVX(t *testing.T) {
	rec := newRecord("read1", "BX:Z:A00C02B03D04")
	valid := HaplotagBX(rec)
	if valid {
		t.Errorf("expected valid=false for barcode containing A00")
	}
	if got, want := string(rec.Desc), "VX:i:0\tBX:Z:A00C02B03D04"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestHaplotagBX_VXOverridesInvalidLookingBarcode(t *testing.T) {
	// VX:i:1 is present and authoritative even though the barcode looks invalid
	rec := newRecord("read1", "VX:i:1\tBX:Z:A00C02B03D04")
	valid := HaplotagBX(rec)
	if !valid {
		t.Errorf("expected valid=true because VX:i:1 is present")
	}
}

func TestHaplotagBX_VXZeroMarksInvalid(t *testing.T) {
	rec := newRecord("read1", "VX:i:0\tBX:Z:A01C02B03D04")
	valid := HaplotagBX(rec)
	if valid {
		t.Errorf("expected valid=false because VX:i:0 overrides a valid-looking barcode")
	}
}

func TestHaplotagBX_NoBXTag(t *testing.T) {
	rec := newRecord("read1", "some other comment")
	valid := HaplotagBX(rec)
	if valid {
		t.Errorf("expected valid=false when no BX tag present")
	}
	want := "some other comment\t" + string(MissingBarcode)
	if got := string(rec.Desc); got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestHaplotagBX_NoBXTagEmptyDesc(t *testing.T) {
	rec := newRecord("read1", "")
	valid := HaplotagBX(rec)
	if valid {
		t.Errorf("expected valid=false when no BX tag present")
	}
	if got, want := string(rec.Desc), string(MissingBarcode); got != want {
		t.Errorf("rec.Desc = %q, want %q (no leading tab expected on empty Desc)", got, want)
	}
}

func TestHaplotagBX_StripsIlluminaOldSuffix(t *testing.T) {
	rec := newRecord("read1/1", "BX:Z:A01C02B03D04")
	HaplotagBX(rec)
	if got, want := string(rec.ID), "read1"; got != want {
		t.Errorf("rec.ID = %q, want %q", got, want)
	}
}

func TestHaplotagBX_StripsIlluminaNewSuffix(t *testing.T) {
	rec := newRecord("read1", "1:N:0:ATCGATCG\tBX:Z:A01C02B03D04")
	HaplotagBX(rec)
	// IlluminaNew's trailing (?:\s|$) consumes the separating tab along with
	// the CASAVA text, so no leftover tab remains before the rebuilt tags.
	if got, want := string(rec.Desc), "VX:i:1\tBX:Z:A01C02B03D04"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestHaplotagBX_RemovesVXAndBXRegardlessOfOrder(t *testing.T) {
	// BX before VX in the description
	rec1 := newRecord("read1", "BX:Z:A01C02B03D04\tVX:i:1\textra")
	valid1 := HaplotagBX(rec1)
	if !valid1 {
		t.Errorf("rec1: expected valid=true")
	}
	if got, want := string(rec1.Desc), "extra\tVX:i:1\tBX:Z:A01C02B03D04"; got != want {
		t.Errorf("rec1.Desc = %q, want %q", got, want)
	}

	// VX before BX: both matches' trailing (?:\s|$) eats the tab that
	// follows them, but the tab preceding "VX" (the one separating it from
	// "extra") is never consumed, so it survives alongside the tab the
	// rebuild step adds — a harmless double tab, not a parsing problem
	// since BX:Z:/VX:i: are matched with \S+ / regexp anchors regardless.
	rec2 := newRecord("read1", "extra\tVX:i:1\tBX:Z:A01C02B03D04")
	valid2 := HaplotagBX(rec2)
	if !valid2 {
		t.Errorf("rec2: expected valid=true")
	}
	if got, want := string(rec2.Desc), "extra\t\tVX:i:1\tBX:Z:A01C02B03D04"; got != want {
		t.Errorf("rec2.Desc = %q, want %q", got, want)
	}
}

func TestTellseq2Std_ValidBarcode(t *testing.T) {
	rec := newRecord("read1:ATCGATCGATCGATCGAT", "")
	valid := Tellseq2Std(rec)
	if !valid {
		t.Errorf("expected valid=true")
	}
	if got, want := string(rec.ID), "read1"; got != want {
		t.Errorf("rec.ID = %q, want %q", got, want)
	}
	if got, want := string(rec.Desc), "VX:i:1\tBX:Z:ATCGATCGATCGATCGAT"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestTellseq2Std_InvalidBarcodeWithN(t *testing.T) {
	rec := newRecord("read1:ATCGNTCGATCGATCGAT", "")
	valid := Tellseq2Std(rec)
	if valid {
		t.Errorf("expected valid=false because barcode contains N")
	}
	if got, want := string(rec.Desc), "VX:i:0\tBX:Z:ATCGNTCGATCGATCGAT"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestTellseq2Std_NoBarcode(t *testing.T) {
	rec := newRecord("read1", "existing comment")
	valid := Tellseq2Std(rec)
	if valid {
		t.Errorf("expected valid=false when no barcode found")
	}
	want := "existing comment\t" + string(MissingBarcode)
	if got := string(rec.Desc); got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
	if got, want := string(rec.ID), "read1"; got != want {
		t.Errorf("rec.ID should be unchanged, got %q, want %q", got, want)
	}
}

func TestStlfr2Std_ValidBarcode(t *testing.T) {
	rec := newRecord("read1#123_456_789", "")
	valid := Stlfr2Std(rec)
	if !valid {
		t.Errorf("expected valid=true")
	}
	if got, want := string(rec.ID), "read1"; got != want {
		t.Errorf("rec.ID = %q, want %q", got, want)
	}
	if got, want := string(rec.Desc), "VX:i:1\tBX:Z:123_456_789"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestStlfr2Std_InvalidZeroComponent(t *testing.T) {
	rec := newRecord("read1#0_456_789", "")
	valid := Stlfr2Std(rec)
	if valid {
		t.Errorf("expected valid=false because barcode has a zero component")
	}
	if got, want := string(rec.Desc), "VX:i:0\tBX:Z:0_456_789"; got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestStlfr2Std_NoBarcode(t *testing.T) {
	rec := newRecord("read1", "")
	valid := Stlfr2Std(rec)
	if valid {
		t.Errorf("expected valid=false when no barcode found")
	}
	if got, want := string(rec.Desc), string(MissingBarcode); got != want {
		t.Errorf("rec.Desc = %q, want %q", got, want)
	}
}

func TestStlfr2Std_StripsIlluminaSuffixes(t *testing.T) {
	// the old-style "/1" suffix only matches at the end of the ID (or before
	// whitespace), so it trails the embedded barcode here.
	rec := newRecord("read1#123_456_789/1", "1:N:0:ATCGATCG")
	valid := Stlfr2Std(rec)
	if !valid {
		t.Errorf("expected valid=true")
	}
	if got, want := string(rec.ID), "read1"; got != want {
		t.Errorf("rec.ID = %q, want %q", got, want)
	}
}
