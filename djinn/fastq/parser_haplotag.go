package fastq

import "github.com/shenwei356/bio/seqio/fastx"

// ----- haplotagging --------------------------------------------------
// Haplotag2Corefq parses a haplotagging-format record — barcode inline as a
// "BX:Z:<bc>" tag in the description — into a CoreFq. ok is false when no
// BX:Z: tag is present
func Haplotag2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	valStart, valEnd, fs, fe, found := findField(desc, BXTAG)
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
	bc := append([]byte(nil), desc[valStart:valEnd]...)
	desc = spliceField(desc, fs, fe)

	return CoreFq{
		ID:       id,
		Seq:      rec.Seq.Seq,
		Qual:     rec.Seq.Qual,
		CASAVA:   casava,
		Barcode:  bc,
		Valid:    !Invalid.Match(bc),
		Comments: desc,
	}, true
}
