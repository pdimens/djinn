package fastq

import (
	"bytes"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func TestToHaplotagging(t *testing.T) {
	id := []byte("read1")
	seq := []byte("ACGT")
	qual := []byte("IIII")
	rec := CoreFq{
		ID:      &id,
		Seq:     &seq,
		Qual:    &qual,
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
	id := []byte("read1")
	seq := []byte("ACGT")
	qual := []byte("IIII")
	rec := CoreFq{
		ID:      &id,
		Seq:     &seq,
		Qual:    &qual,
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
	id := []byte("read1")
	seq := []byte("ACGT")
	qual := []byte("IIII")
	rec := CoreFq{
		ID:      &id,
		Seq:     &seq,
		Qual:    &qual,
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
	id := []byte("read1")
	seq := []byte("ACGT")
	qual := []byte("IIII")
	rec := CoreFq{
		ID:      &id,
		Seq:     &seq,
		Qual:    &qual,
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
	id := []byte("read1")
	seq := []byte("ACGTACGT")
	qual := []byte("IIIIIIII")
	rec := CoreFq{
		ID:      &id,
		Seq:     &seq,
		Qual:    &qual,
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
	id := []byte("read1")
	seq := []byte("ACGTACGT")
	qual := []byte("IIIIIIII")
	rec := CoreFq{
		ID:      &id,
		Seq:     &seq,
		Qual:    &qual,
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

// TestHaplotag2Corefq_KnownLimitations documents the current, incomplete
// behavior of Haplotag2Corefq: (1) it hardcodes Barcode to "AT" regardless
// of the input record — it doesn't actually parse a barcode out of rec —
// and (2) it slices rec.Desc[:2] with no length check, which panics for any
// record whose Desc is shorter than 2 bytes. It also unconditionally
// dereferences rec.Seq, so a Record without a populated Seq (as produced by
// some code paths, or in a minimal test fixture) panics too. This function
// is not wired up to any caller yet elsewhere in the repo; these tests
// exist to pin down and flag its current (buggy/stub) behavior rather than
// to validate it as correct.
func TestHaplotag2Corefq_KnownLimitations(t *testing.T) {
	desc := []byte("XY some description")
	rec := &fastx.Record{
		ID:   []byte("read1"),
		Desc: desc,
		Seq:  &seq.Seq{Seq: []byte("ACGT"), Qual: []byte("IIII")},
	}
	core := Haplotag2Corefq(rec)

	if got, want := string(core.CASAVA), "XY"; got != want {
		t.Errorf("CASAVA = %q, want %q (first 2 bytes of Desc)", got, want)
	}
	// This is the bug: Barcode is always "AT" no matter what the record
	// actually contains.
	if got, want := string(core.Barcode), "AT"; got != want {
		t.Errorf("Barcode = %q, want %q (hardcoded stub value)", got, want)
	}

	t.Run("panics on short Desc", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Errorf("expected Haplotag2Corefq to panic on a Desc shorter than 2 bytes (known bug, not fixed)")
			}
		}()
		shortRec := &fastx.Record{
			ID:   []byte("read1"),
			Desc: []byte("X"),
			Seq:  &seq.Seq{Seq: []byte("ACGT"), Qual: []byte("IIII")},
		}
		_ = Haplotag2Corefq(shortRec)
	})

	t.Run("panics on nil Seq", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Errorf("expected Haplotag2Corefq to panic when rec.Seq is nil (known bug, not fixed)")
			}
		}()
		nilSeqRec := &fastx.Record{ID: []byte("read1"), Desc: desc}
		_ = Haplotag2Corefq(nilSeqRec)
	})
}
