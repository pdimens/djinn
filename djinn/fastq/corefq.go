package fastq

import (
	"bytes"

	"github.com/shenwei356/bio/seqio/fastx"
)

type CoreFq struct {
	ID       *[]byte // the ID/Name of a fastq record
	Qual     *[]byte // the PHRED quality score of the fastq record
	Seq      *[]byte // the sequence of a fastq record
	Comments *[]byte // the comments of a fastq record header
	CASAVA   []byte  // the CASAVA part of a read ID
	Barcode  []byte  // the isolated barcode of a record
}

// Process a fastx.Record into a CoreFq type
func Haplotag2Corefq(rec *fastx.Record) CoreFq {

	return CoreFq{
		ID:       &rec.ID,
		Qual:     &rec.Seq.Qual,
		Seq:      &rec.Seq.Seq,
		Comments: &rec.Desc,
		CASAVA:   rec.Desc[:2],
		Barcode:  []byte{'A', 'T'},
	}
}

// ---- Conversions -------------------------------------------------

// Convert the fastq record to haplotagging format and write it the write buffer
func ToHaplotagging(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(*rec.ID)
	writer.WriteByte('/')
	// only use the first byte of CASAVA, which will either be 1 or 2
	writer.WriteByte(rec.CASAVA[0])
	writer.WriteByte(TabSep)
	//TODO ONLY WRITE SAM-COMPLIANT COMMENTS
	writer.Write([]byte{'B', 'X', ':', 'Z', ':'})
	writer.Write(rec.Barcode)
	writer.WriteByte(Newline)
	writer.Write(*rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(*rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to tellseq format and write it the write buffer
func ToTellseq(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(*rec.ID)
	writer.WriteByte(':')
	writer.Write(rec.Barcode)
	writer.WriteByte(TabSep)
	// only use the first byte of CASAVA, which will either be 1 or 2
	if n := len(rec.CASAVA); n == 1 {
		writer.Write([]byte{rec.CASAVA[0], ':', 'N', ':', 'A', 'T', 'C', 'G'})
	} else {
		writer.Write(rec.CASAVA)
	}
	writer.WriteByte(Newline)
	writer.Write(*rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(*rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to haplotagging format and write it the write buffer
func ToStlfr(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(*rec.ID)
	writer.WriteByte('#')
	writer.Write(rec.Barcode)
	writer.WriteByte(TabSep)
	// only use the first byte of CASAVA, which will either be 1 or 2
	if n := len(rec.CASAVA); n == 1 {
		writer.Write([]byte{rec.CASAVA[0], ':', 'N', ':', 'A', 'T', 'C', 'G'})
	} else {
		writer.Write(rec.CASAVA)
	}
	writer.WriteByte(Newline)
	writer.Write(*rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(*rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to 10X format and write it the write buffer
func ToTenX(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(*rec.ID)
	writer.WriteByte('#')
	writer.Write(rec.Barcode)
	writer.WriteByte(TabSep)
	// only use the first byte of CASAVA, which will either be 1 or 2
	if n := len(rec.CASAVA); n == 1 {
		writer.Write([]byte{rec.CASAVA[0], ':', 'N', ':', 'A', 'T', 'C', 'G'})
	} else {
		writer.Write(rec.CASAVA)
	}
	writer.WriteByte(Newline)
	// only write barcode as R1 seq prefix
	if rec.CASAVA[0] == '1' {
		writer.Write(rec.Barcode)
		writer.Write(*rec.Seq)
		writer.WriteByte(Newline)
		writer.WriteByte(PlusSign)
		writer.WriteByte(Newline)
		for range rec.Barcode {
			writer.WriteByte('I')
		}
		writer.Write(*rec.Qual)
	} else {
		writer.Write(*rec.Seq)
		writer.WriteByte(Newline)
		writer.WriteByte(PlusSign)
		writer.WriteByte(Newline)
		writer.Write(*rec.Qual)
	}
	writer.WriteByte(Newline)
}
