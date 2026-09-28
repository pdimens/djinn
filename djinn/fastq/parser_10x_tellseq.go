package fastq

import (
	"bytes"

	"github.com/shenwei356/bio/seqio/fastx"
)

func isValidNuc(seq []byte) bool {
	if bytes.IndexByte(seq, 'N') != -1 {
		return false // Found 'N'
	}
	return true
}

func TenX2Corefq(rec *fastx.Record, bc_map map[string]struct{}) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)
	bc := rec.Seq.Seq[0:16]
	if _, ok := bc_map[string(bc)]; ok {
		bc = append([]byte(nil), rec.Seq.Seq[0:16]...)
		rec.Seq.Seq = append([]byte(nil), rec.Seq.Seq[16:]...)
		rec.Seq.Qual = append([]byte(nil), rec.Seq.Qual[16:]...)
		return CoreFq{
			ID:       id,
			Seq:      rec.Seq.Seq,
			Qual:     rec.Seq.Qual,
			CASAVA:   casava,
			Barcode:  bc,
			Valid:    isValidNuc(bc),
			Comments: desc,
		}, true
	}
	return CoreFq{
		ID:       id,
		Seq:      rec.Seq.Seq,
		Qual:     rec.Seq.Qual,
		CASAVA:   casava,
		Barcode:  nil,
		Valid:    false,
		Comments: desc,
	}, false
}

// ---- tellseq -----------------------------

// findTellseqBarcode locates a ":<ACGTN+>" suffix at the end of id
// (Tellseq's inline barcode format), mirroring Tellseq.FindSubmatchIndex
// applied to id (which, containing no whitespace, can only match at the
// very end). Returns the barcode's [start,end) and the index of the
// leading ':' to splice from.
func findTellseqBarcode(id []byte) (start, end, colon int, ok bool) {
	n := len(id)
	i := n
	for i > 0 && isBase(id[i-1]) {
		i--
	}
	if i == n || i == 0 || id[i-1] != ':' {
		return 0, 0, 0, false
	}
	return i, n, i - 1, true
}

// Tellseq2Corefq parses a TELLseq-format record — barcode inline at the end
// of the ID as ":<bases>" — into a CoreFq. ok is false when no such suffix
// is present
func Tellseq2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	start, end, colon, found := findTellseqBarcode(id)
	if !found {
		return CoreFq{
			ID:       id,
			Seq:      rec.Seq.Seq,
			Qual:     rec.Seq.Qual,
			CASAVA:   casava,
			Barcode:  nil,
			Valid:    false,
			Comments: desc,
		}, false
	}
	bc := append([]byte(nil), id[start:end]...)
	id = append(id[:colon], id[end:]...)

	return CoreFq{
		ID:       id,
		Seq:      rec.Seq.Seq,
		Qual:     rec.Seq.Qual,
		CASAVA:   casava,
		Barcode:  bc,
		Valid:    isValidNuc(bc),
		Comments: desc,
	}, true
}
