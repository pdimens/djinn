package fastq

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"math"
	"strconv"

	"github.com/shenwei356/bio/seqio/fastx"
)

type CoreFq struct {
	ID       *[]byte // the ID/Name of a fastq record
	Qual     *[]byte // the PHRED quality score of the fastq record
	Seq      *[]byte // the sequence of a fastq record
	Comments []byte  // the indices of the comments of a fastq record header
	CASAVA   []byte  // the CASAVA part of a read ID
	Barcode  []byte  // the isolated barcode of a record
}

// ----- haplotagging ----------------------------------

// bxTag returns the "BX:Z:..." field from an aux byte slice, up to the next
// tab or end of slice. Returns nil if not present.
func bxTagIndices(aux []byte) [2]int {
	var indices [2]int
	i := bytes.Index(aux, []byte{'B', 'X', ':', 'Z', ':'})
	if i == -1 {
		return indices
	}
	j := bytes.Index(aux[i:], []byte{'\t'})
	if j != -1 {
		indices[1] = len(aux)
	} else {
		indices[1] = j
	}
	return indices
}

// Process a haplotagging fastx.Record into a CoreFq type
func Haplotag2Corefq(rec *fastx.Record) CoreFq {
	return CoreFq{
		ID:       &rec.ID,
		Qual:     &rec.Seq.Qual,
		Seq:      &rec.Seq.Seq,
		Comments: rec.Desc,
		CASAVA:   rec.Desc[len(rec.Desc)-2:],
		Barcode:  bxTag(rec.Desc),
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

// --- write directly to BAM without sam.Record --------
var baseCode = func() (t [256]byte) {
	for i := range t {
		t[i] = 15
	} // N
	for k, v := range map[byte]byte{'=': 0, 'A': 1, 'C': 2, 'M': 3, 'G': 4, 'R': 5, 'S': 6, 'V': 7,
		'T': 8, 'W': 9, 'Y': 10, 'H': 11, 'K': 12, 'D': 13, 'B': 14, 'N': 15} {
		t[k] = v
	}
	return
}()

func reg2bin(beg, end int32) uint16 {
	end--
	switch {
	case beg>>14 == end>>14:
		return uint16(((1<<15)-1)/7 + (beg >> 14))
	case beg>>17 == end>>17:
		return uint16(((1<<12)-1)/7 + (beg >> 17))
	case beg>>20 == end>>20:
		return uint16(((1<<9)-1)/7 + (beg >> 20))
	case beg>>23 == end>>23:
		return uint16(((1<<6)-1)/7 + (beg >> 23))
	case beg>>26 == end>>26:
		return uint16(((1<<3)-1)/7 + (beg >> 26))
	}
	return 0
}

// encodeAux converts SAM-text aux fields ("XX:i:123\tYY:Z:foo ...", tab or space
// sep) into BAM binary aux bytes. Numeric types collapsed to 'i'(int32) and
// 'f'(float32) for simplicity — not size-optimal like htslib, but correct.
// B/H array tags: unsupported, returns error.
func encodeAux(comments []byte) ([]byte, error) {
	if len(comments) == 0 {
		return nil, nil
	}
	buf := make([]byte, 0, len(comments)+8)
	fields := bytes.FieldsFunc(comments, func(r rune) bool { return r == '\t' || r == ' ' })
	for _, f := range fields {
		if len(f) < 5 || f[2] != ':' || f[4] != ':' {
			return nil, fmt.Errorf("bad aux field: %q", f)
		}
		tag, typ, val := f[0:2], f[3], f[5:]
		buf = append(buf, tag...)
		switch typ {
		case 'A':
			buf = append(buf, 'A', val[0])
		case 'i':
			n, err := strconv.ParseInt(string(val), 10, 32)
			if err != nil {
				return nil, err
			}
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], uint32(int32(n)))
			buf = append(buf, 'i')
			buf = append(buf, b[:]...)
		case 'f':
			fl, err := strconv.ParseFloat(string(val), 32)
			if err != nil {
				return nil, err
			}
			var b [4]byte
			binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(fl)))
			buf = append(buf, 'f')
			buf = append(buf, b[:]...)
		case 'Z':
			buf = append(buf, 'Z')
			buf = append(buf, val...)
			buf = append(buf, 0)
		case 'H':
			return nil, fmt.Errorf("H aux type unsupported: %q", tag)
		case 'B':
			return nil, fmt.Errorf("B aux type unsupported: %q", tag)
		default:
			return nil, fmt.Errorf("unknown aux type %c in %q", typ, tag)
		}
	}
	return buf, nil
}

var negOne int32 = -1

// Convert to and write to buffer as unmapped BAM record, direct bytes, no sam.Record intermediate.
func (c *CoreFq) WriteBam(w *bufio.Writer) error {
	id, seq, qual := *c.ID, *c.Seq, *c.Qual
	lReadName := len(id) + 1
	lSeq := len(seq)

	aux, err := encodeAux(*c.Comments)
	if err != nil {
		return err
	}

	var flag uint16 = 4 // unmapped
	if len(c.CASAVA) > 0 {
		flag |= 1 // paired
		switch c.CASAVA[0] {
		case '1':
			flag |= 0x40 // first in pair
		case '2':
			flag |= 0x80 // last in pair
		}
	}

	blockSize := 32 + lReadName + (lSeq+1)/2 + lSeq + len(aux)

	var hdr [36]byte
	binary.LittleEndian.PutUint32(hdr[0:4], uint32(blockSize))
	binary.LittleEndian.PutUint32(hdr[4:8], uint32(negOne))
	binary.LittleEndian.PutUint32(hdr[8:12], uint32(negOne))
	hdr[12] = byte(lReadName)
	hdr[13] = 0
	binary.LittleEndian.PutUint16(hdr[14:16], reg2bin(-1, 0))
	binary.LittleEndian.PutUint16(hdr[16:18], 0)
	binary.LittleEndian.PutUint16(hdr[18:20], flag)
	binary.LittleEndian.PutUint32(hdr[20:24], uint32(lSeq))
	binary.LittleEndian.PutUint32(hdr[24:28], uint32(negOne))
	binary.LittleEndian.PutUint32(hdr[28:32], uint32(negOne))
	binary.LittleEndian.PutUint32(hdr[32:36], 0)

	if _, err := w.Write(hdr[:]); err != nil {
		return err
	}
	if _, err := w.Write(id); err != nil {
		return err
	}
	if err := w.WriteByte(0); err != nil {
		return err
	}

	for i := 0; i < lSeq; i += 2 {
		hi := baseCode[seq[i]]
		var lo byte
		if i+1 < lSeq {
			lo = baseCode[seq[i+1]]
		}
		if err := w.WriteByte(hi<<4 | lo); err != nil {
			return err
		}
	}

	if len(qual) == lSeq {
		for i := range qual {
			if err := w.WriteByte(qual[i] - 33); err != nil {
				return err
			}
		}
	} else {
		for range lSeq {
			if err := w.WriteByte(0xff); err != nil {
				return err
			}
		}
	}

	if len(aux) > 0 {
		if _, err := w.Write(aux); err != nil {
			return err
		}
	}
	return nil
}
