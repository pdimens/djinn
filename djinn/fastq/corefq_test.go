package fastq

import (
	"bytes"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func TestToHaplotagging(t *testing.T) {
	rec := CoreFq{
		ID:      []byte("read1"),
		Seq:     []byte("ACGT"),
		Qual:    []byte("IIII"),
		CASAVA:  []byte("1"),
		Barcode: []byte("A01C02B03D04"),
	}
	var buf bytes.Buffer
	ToHaplotagging(&rec, &buf)

	want := "@read1/1\tBX:Z:A01C02B03D04\nACGT\n+\nIIII\n"
	if got := buf.String(); got != want {
		t.Errorf("ToHaplotagging output = %q, want %q", got, want)
	}
}

func TestToTellseq_SingleByteCASAVA(t *testing.T) {
	rec := CoreFq{
		ID:      []byte("read1"),
		Seq:     []byte("ACGT"),
		Qual:    []byte("IIII"),
		CASAVA:  []byte("2"),
		Barcode: []byte("ATCGATCG"),
	}
	var buf bytes.Buffer
	ToTellseq(&rec, &buf)

	want := "@read1:ATCGATCG\t2:N:ATCG\nACGT\n+\nIIII\n"
	if got := buf.String(); got != want {
		t.Errorf("ToTellseq output = %q, want %q", got, want)
	}
}

func TestToTellseq_FullCASAVA(t *testing.T) {
	rec := CoreFq{
		ID:      []byte("read1"),
		Seq:     []byte("ACGT"),
		Qual:    []byte("IIII"),
		CASAVA:  []byte("1:N:0:ATAG"),
		Barcode: []byte("ATCGATCG"),
	}
	var buf bytes.Buffer
	ToTellseq(&rec, &buf)

	want := "@read1:ATCGATCG\t1:N:0:ATAG\nACGT\n+\nIIII\n"
	if got := buf.String(); got != want {
		t.Errorf("ToTellseq output = %q, want %q", got, want)
	}
}

func TestToStlfr(t *testing.T) {
	rec := CoreFq{
		ID:      []byte("read1"),
		Seq:     []byte("ACGT"),
		Qual:    []byte("IIII"),
		CASAVA:  []byte("1"),
		Barcode: []byte("123_456_789"),
	}
	var buf bytes.Buffer
	ToStlfr(&rec, &buf)

	want := "@read1#123_456_789\t1:N:ATCG\nACGT\n+\nIIII\n"
	if got := buf.String(); got != want {
		t.Errorf("ToStlfr output = %q, want %q", got, want)
	}
}

func TestToTenX_R1PrependsBarcodeAndFillerQual(t *testing.T) {
	rec := CoreFq{
		ID:      []byte("read1"),
		Seq:     []byte("ACGTACGT"),
		Qual:    []byte("IIIIIIII"),
		CASAVA:  []byte("1"),
		Barcode: []byte("ATCGATCGATCGATCG"), // 16bp 10x barcode
	}
	var buf bytes.Buffer
	ToTenX(&rec, &buf)

	want := "@read1#ATCGATCGATCGATCG\t1:N:ATCG\n" +
		"ATCGATCGATCGATCGACGTACGT\n+\n" +
		"IIIIIIIIIIIIIIIIIIIIIIII\n" // 16 'I' filler + original 8 qual chars
	if got := buf.String(); got != want {
		t.Errorf("ToTenX (R1) output = %q, want %q", got, want)
	}
}

func TestToTenX_R2LeavesSeqUnmodified(t *testing.T) {
	rec := CoreFq{
		ID:      []byte("read1"),
		Seq:     []byte("ACGTACGT"),
		Qual:    []byte("IIIIIIII"),
		CASAVA:  []byte("2"),
		Barcode: []byte("ATCGATCGATCGATCG"),
	}
	var buf bytes.Buffer
	ToTenX(&rec, &buf)

	want := "@read1#ATCGATCGATCGATCG\t2:N:ATCG\nACGTACGT\n+\nIIIIIIII\n"
	if got := buf.String(); got != want {
		t.Errorf("ToTenX (R2) output = %q, want %q", got, want)
	}
}

// ---- constructors --------------------------------------------------

func mkRecord(id, desc, seq_, qual string) *fastx.Record {
	return &fastx.Record{
		ID:   []byte(id),
		Desc: []byte(desc),
		Seq:  &seq.Seq{Seq: []byte(seq_), Qual: []byte(qual)},
	}
}

func TestHaplotag2Corefq(t *testing.T) {
	t.Run("BX and new-style CASAVA", func(t *testing.T) {
		rec := mkRecord("read1", "1:N:0:ATCG BX:Z:A01C02B03D04", "ACGT", "IIII")
		core, ok := Haplotag2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if got, want := string(core.ID), "read1"; got != want {
			t.Errorf("ID = %q, want %q", got, want)
		}
		if got, want := string(core.CASAVA), "1:N:0:ATCG"; got != want {
			t.Errorf("CASAVA = %q, want %q", got, want)
		}
		if got, want := string(core.Barcode), "A01C02B03D04"; got != want {
			t.Errorf("Barcode = %q, want %q", got, want)
		}
		if !core.Valid {
			t.Error("expected Valid=true for a clean barcode")
		}
		if got, want := string(core.Comments), ""; got != want {
			t.Errorf("Comments = %q, want %q (BX and CASAVA both removed)", got, want)
		}
		if got, want := string(core.Seq), "ACGT"; got != want {
			t.Errorf("Seq = %q, want %q", got, want)
		}
	})

	t.Run("old-style CASAVA suffix on ID", func(t *testing.T) {
		rec := mkRecord("read1/2", "BX:Z:A01C02B03D04", "ACGT", "IIII")
		core, ok := Haplotag2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if got, want := string(core.ID), "read1"; got != want {
			t.Errorf("ID = %q, want %q", got, want)
		}
		if got, want := string(core.CASAVA), "2"; got != want {
			t.Errorf("CASAVA = %q, want %q", got, want)
		}
	})

	t.Run("invalid barcode (contains N)", func(t *testing.T) {
		rec := mkRecord("read1", "BX:Z:A01C02BNND04", "ACGT", "IIII")
		core, ok := Haplotag2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if core.Valid {
			t.Error("expected Valid=false for a barcode containing N")
		}
	})

	t.Run("other aux fields survive in Comments", func(t *testing.T) {
		rec := mkRecord("read1", "XX:i:5 BX:Z:A01C02B03D04 YY:Z:foo", "ACGT", "IIII")
		core, ok := Haplotag2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if got, want := string(core.Comments), "XX:i:5 YY:Z:foo"; got != want {
			t.Errorf("Comments = %q, want %q", got, want)
		}
	})

	t.Run("no BX tag", func(t *testing.T) {
		rec := mkRecord("read1", "some other text", "ACGT", "IIII")
		_, ok := Haplotag2Corefq(rec)
		if ok {
			t.Error("expected ok=false when no BX:Z: tag is present")
		}
	})
}

func TestTellseq2Corefq(t *testing.T) {
	t.Run("barcode and old-style CASAVA", func(t *testing.T) {
		rec := mkRecord("read1:ATCGATCG/1", "", "ACGT", "IIII")
		core, ok := Tellseq2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if got, want := string(core.ID), "read1"; got != want {
			t.Errorf("ID = %q, want %q", got, want)
		}
		if got, want := string(core.Barcode), "ATCGATCG"; got != want {
			t.Errorf("Barcode = %q, want %q", got, want)
		}
		if got, want := string(core.CASAVA), "1"; got != want {
			t.Errorf("CASAVA = %q, want %q", got, want)
		}
		if !core.Valid {
			t.Error("expected Valid=true")
		}
	})

	t.Run("no barcode", func(t *testing.T) {
		rec := mkRecord("read1", "", "ACGT", "IIII")
		_, ok := Tellseq2Corefq(rec)
		if ok {
			t.Error("expected ok=false when ID has no :ACGTN+ suffix")
		}
	})

	t.Run("invalid barcode", func(t *testing.T) {
		rec := mkRecord("read1:ATCGNNCG", "", "ACGT", "IIII")
		core, ok := Tellseq2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if core.Valid {
			t.Error("expected Valid=false for a barcode containing N")
		}
	})
}

func TestStlfr2Corefq(t *testing.T) {
	t.Run("barcode and new-style CASAVA", func(t *testing.T) {
		rec := mkRecord("read1#123_456_789", "2:N:0:ATCG", "ACGT", "IIII")
		core, ok := Stlfr2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if got, want := string(core.ID), "read1"; got != want {
			t.Errorf("ID = %q, want %q", got, want)
		}
		if got, want := string(core.Barcode), "123_456_789"; got != want {
			t.Errorf("Barcode = %q, want %q", got, want)
		}
		if got, want := string(core.CASAVA), "2:N:0:ATCG"; got != want {
			t.Errorf("CASAVA = %q, want %q", got, want)
		}
		if !core.Valid {
			t.Error("expected Valid=true")
		}
	})

	t.Run("invalid barcode (zero index)", func(t *testing.T) {
		rec := mkRecord("read1#123_0_789", "", "ACGT", "IIII")
		core, ok := Stlfr2Corefq(rec)
		if !ok {
			t.Fatal("expected ok=true")
		}
		if core.Valid {
			t.Error("expected Valid=false for a barcode with a zero index")
		}
	})

	t.Run("no barcode", func(t *testing.T) {
		rec := mkRecord("read1", "", "ACGT", "IIII")
		_, ok := Stlfr2Corefq(rec)
		if ok {
			t.Error("expected ok=false when ID has no #n_n_n suffix")
		}
	})
}

// aliasing regression: constructors must not corrupt Desc/ID by reusing
// their backing array across the splice of an earlier-extracted substring.
func TestCorefqConstructors_NoAliasingCorruption(t *testing.T) {
	rec := mkRecord("readAAAA:ATCGATCGATCGATCG", "1:N:0:ATCGATCGATCG", "ACGT", "IIII")
	core, ok := Tellseq2Corefq(rec)
	if !ok {
		t.Fatal("expected ok=true")
	}
	// force a reallocation-triggering append to make sure Barcode/CASAVA
	// were independent copies, not still-aliased views into rec.ID/rec.Desc
	rec.ID = append(rec.ID, 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X')
	if got, want := string(core.Barcode), "ATCGATCGATCGATCG"; got != want {
		t.Errorf("Barcode corrupted after mutating rec.ID: got %q, want %q", got, want)
	}
	if got, want := string(core.CASAVA), "1:N:0:ATCGATCGATCG"; got != want {
		t.Errorf("CASAVA corrupted: got %q, want %q", got, want)
	}
}
