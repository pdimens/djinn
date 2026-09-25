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

// CoreFq is a zero-copy, per-record view over the parts of a FASTQ record
// that matter for linked-read barcode handling: sequence/quality, the
// paired-end (CASAVA) marker, the isolated barcode, and (for direct BAM
// output) whatever SAM-style aux text is left in the description once the
// barcode and CASAVA marker have been pulled out of it.
//
// Every []byte field aliases memory owned by the *fastx.Record (or, more
// precisely, the fastx.Reader that produced it — see below) CoreFq was
// built from; building one does not allocate beyond copying out the short
// Barcode/CASAVA substrings that must survive being spliced out of ID/Desc.
//
// Lifetime: fastx.Reader.Read() reuses a single Record and its internal
// buffers across calls — it does not allocate a fresh one each time. A
// CoreFq is therefore only valid until the next call to that reader's
// Read(): build it, use it (convert/write it), and let it go before
// reading the next record. Never hand a CoreFq, or any of its fields, to
// another goroutine while the reader keeps advancing — that reintroduces
// the same class of data corruption bug this design is built to avoid (see
// HaplotagBX's bc-aliasing fix in processfq.go). For real concurrency, fan
// out over fastx.Reader.ChunkChan(bufferSize, chunkSize) instead: it calls
// Record.Clone() before a record crosses the channel, so each chunk owns
// independent memory and a CoreFq built from it inside a worker goroutine
// is safe on its own.
type CoreFq struct {
	ID       []byte // read ID, CASAVA suffix removed
	Seq      []byte // nucleotide sequence
	Qual     []byte // PHRED quality string
	CASAVA   []byte // paired-end marker: "1"/"2" (old style) or "1:N:0:ATCG" (new style); nil if none was found
	Barcode  []byte // isolated linked-read barcode; nil if none was found
	Valid    bool   // whether Barcode passes the Invalid pattern; meaningless when Barcode is nil
	Comments []byte // remaining description text (SAM-style aux tokens) once CASAVA/Barcode are removed
}

// ----- shared scanning helpers -------------------------------------------
//
// These replace regexp.FindSubmatchIndex-based matching (IlluminaNew,
// IlluminaOld, StdBx, Tellseq, Stlfr in barcodes.go) with single-pass byte
// scans that match the exact same rules. Regexp is fine at the call volumes
// those functions see elsewhere, but CoreFq is meant to be built once per
// record at tens-of-millions-of-records scale, so it avoids the engine
// overhead entirely.

func isSpace(b byte) bool { return b == ' ' || b == '\t' }

func isBase(b byte) bool {
	switch b {
	case 'A', 'C', 'G', 'T', 'N':
		return true
	}
	return false
}

func isDigit(b byte) bool { return b >= '0' && b <= '9' }

func isAlnum(b byte) bool {
	return b >= 'A' && b <= 'Z' || b >= 'a' && b <= 'z' || isDigit(b)
}

// nextField returns the [start,end) range of the next whitespace-delimited
// field in data at or after idx, and the offset to resume scanning from.
// ok is false once no more fields remain.
func nextField(data []byte, idx int) (start, end, next int, ok bool) {
	n := len(data)
	for idx < n && isSpace(data[idx]) {
		idx++
	}
	if idx >= n {
		return 0, 0, n, false
	}
	start = idx
	for idx < n && !isSpace(data[idx]) {
		idx++
	}
	return start, idx, idx, true
}

// findField scans desc field-by-field for one that starts with prefix, and
// returns the [valStart,valEnd) range of the value after prefix along with
// the [fieldStart,fieldEnd) range of the whole field (for splicing it back
// out of desc). ok is false if no field has that prefix.
func findField(desc, prefix []byte) (valStart, valEnd, fieldStart, fieldEnd int, ok bool) {
	for i := 0; ; {
		fs, fe, next, more := nextField(desc, i)
		if !more {
			return 0, 0, 0, 0, false
		}
		i = next
		if fe-fs > len(prefix) && bytes.HasPrefix(desc[fs:fe], prefix) {
			return fs + len(prefix), fe, fs, fe, true
		}
	}
}

// spliceField removes data[start:end) in place, absorbing one adjacent
// separator byte if present so a stray double-separator isn't left behind.
func spliceField(data []byte, start, end int) []byte {
	if end < len(data) && isSpace(data[end]) {
		end++
	} else if start > 0 && isSpace(data[start-1]) {
		start--
	}
	return append(data[:start], data[end:]...)
}

// isIlluminaNewField reports whether f matches the new-style CASAVA marker
// "[12]:[YN]:<digits>:<alnum+>" (e.g. "1:N:0:ATCG"), mirroring IlluminaNew.
func isIlluminaNewField(f []byte) bool {
	if len(f) < 7 || (f[0] != '1' && f[0] != '2') || f[1] != ':' {
		return false
	}
	if f[2] != 'Y' && f[2] != 'N' {
		return false
	}
	if f[3] != ':' {
		return false
	}
	i := 4
	digitsStart := i
	for i < len(f) && isDigit(f[i]) {
		i++
	}
	if i == digitsStart || i >= len(f) || f[i] != ':' {
		return false
	}
	i++
	if i >= len(f) {
		return false
	}
	for ; i < len(f); i++ {
		if !isAlnum(f[i]) {
			return false
		}
	}
	return true
}

// findIlluminaNew locates a new-style CASAVA field in desc, mirroring
// IlluminaNew.FindIndex(desc).
func findIlluminaNew(desc []byte) (start, end int, ok bool) {
	for i := 0; ; {
		fs, fe, next, more := nextField(desc, i)
		if !more {
			return 0, 0, false
		}
		i = next
		if isIlluminaNewField(desc[fs:fe]) {
			return fs, fe, true
		}
	}
}

// findIlluminaOld locates an old-style "/1" or "/2" CASAVA suffix at the
// end of id, mirroring IlluminaOld.FindIndex(id) — id never contains
// whitespace (fastx.Reader splits it from desc at the first whitespace),
// so the only place IlluminaOld's "/[12](?:\s|$)" can match is the very end.
func findIlluminaOld(id []byte) (start, end int, ok bool) {
	n := len(id)
	if n < 2 || id[n-2] != '/' {
		return 0, 0, false
	}
	if id[n-1] != '1' && id[n-1] != '2' {
		return 0, 0, false
	}
	return n - 2, n, true
}

// splitCASAVA removes a paired-end/CASAVA marker from id or desc — new
// style lives in desc, old style at the end of id — and returns the
// (possibly shortened) id/desc plus a copy of the extracted marker. The
// copy is necessary: it's read back out after the in-place splice below,
// which reuses the same backing array (the same hazard HaplotagBX's bc
// copy guards against in processfq.go).
func splitCASAVA(id, desc []byte) (newID, newDesc, casava []byte) {
	if s, e, ok := findIlluminaNew(desc); ok {
		casava = append([]byte(nil), desc[s:e]...)
		return id, spliceField(desc, s, e), casava
	}
	if s, e, ok := findIlluminaOld(id); ok {
		casava = append([]byte(nil), id[s+1:e]...) // digit only, "/" dropped
		return append(id[:s], id[e:]...), desc, casava
	}
	return id, desc, nil
}

// ----- haplotagging --------------------------------------------------

// Haplotag2Corefq parses a haplotagging-format record — barcode inline as a
// "BX:Z:<bc>" tag in the description — into a CoreFq. ok is false when no
// BX:Z: tag is present, in which case the returned CoreFq is the zero
// value.
func Haplotag2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	valStart, valEnd, fs, fe, found := findField(desc, BXTAG)
	if !found {
		return CoreFq{}, false
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

// ----- tellseq ---------------------------------------------------------

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
// is present, in which case the returned CoreFq is the zero value.
func Tellseq2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	start, end, colon, found := findTellseqBarcode(id)
	if !found {
		return CoreFq{}, false
	}
	bc := append([]byte(nil), id[start:end]...)
	id = append(id[:colon], id[end:]...)

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

// ----- stlfr -------------------------------------------------------------

// isStlfrBody reports whether b matches "[0-9]+_[0-9]+_[0-9]+" exactly.
func isStlfrBody(b []byte) bool {
	parts, start, digits := 0, 0, false
	for i := 0; i <= len(b); i++ {
		if i == len(b) || b[i] == '_' {
			if !digits || i == start {
				return false
			}
			parts++
			start = i + 1
			digits = false
			continue
		}
		if !isDigit(b[i]) {
			return false
		}
		digits = true
	}
	return parts == 3
}

// findStlfrBarcode locates a "#<n>_<n>_<n>" suffix in id, mirroring
// Stlfr.FindSubmatchIndex applied to id.
func findStlfrBarcode(id []byte) (start, end, hash int, ok bool) {
	for i := 0; i < len(id); i++ {
		if id[i] != '#' {
			continue
		}
		if body := id[i+1:]; isStlfrBody(body) {
			return i + 1, len(id), i, true
		}
	}
	return 0, 0, 0, false
}

// Stlfr2Corefq parses an stLFR-format record — barcode inline at the end of
// the ID as "#n_n_n" — into a CoreFq. ok is false when no such suffix is
// present, in which case the returned CoreFq is the zero value.
func Stlfr2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	start, end, hash, found := findStlfrBarcode(id)
	if !found {
		return CoreFq{}, false
	}
	bc := append([]byte(nil), id[start:end]...)
	id = append(id[:hash], id[end:]...)

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

// ---- Conversions -------------------------------------------------
//
// These assume rec.CASAVA is non-empty (as every constructor above
// populates it when a marker is present); a record with no CASAVA marker
// at all is a pre-existing gap in this format, not one introduced here.

// Convert the fastq record to haplotagging format and write it the write buffer
func ToHaplotagging(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
	writer.WriteByte('/')
	// only use the first byte of CASAVA, which will either be 1 or 2
	writer.WriteByte(rec.CASAVA[0])
	writer.WriteByte(TabSep)
	//TODO ONLY WRITE SAM-COMPLIANT COMMENTS
	writer.Write([]byte{'B', 'X', ':', 'Z', ':'})
	writer.Write(rec.Barcode)
	writer.WriteByte(Newline)
	writer.Write(rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to tellseq format and write it the write buffer
func ToTellseq(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
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
	writer.Write(rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to haplotagging format and write it the write buffer
func ToStlfr(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
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
	writer.Write(rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to 10X format and write it the write buffer
func ToTenX(rec *CoreFq, writer *bytes.Buffer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
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
		writer.Write(rec.Seq)
		writer.WriteByte(Newline)
		writer.WriteByte(PlusSign)
		writer.WriteByte(Newline)
		for range rec.Barcode {
			writer.WriteByte('I')
		}
		writer.Write(rec.Qual)
	} else {
		writer.Write(rec.Seq)
		writer.WriteByte(Newline)
		writer.WriteByte(PlusSign)
		writer.WriteByte(Newline)
		writer.Write(rec.Qual)
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
	id, seq, qual := c.ID, c.Seq, c.Qual
	lReadName := len(id) + 1
	lSeq := len(seq)

	aux, err := encodeAux(c.Comments)
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
