package xam

// preproc_barcodes - renames and records demultiplexed linked-read barcodes
// coming out of Pheniqs. Reads a Pheniqs JSON config and a BAM file, rewrites
// BX and VX tags on every record that has an RX tag, and writes R1/R2 FASTQ
// output files compressed with pgzip.
//
// Usage:
//
//	preproc_barcodes [--threads N] <pheniqs.json> <input.bam> <out_R1.fq.gz> <out_R2.fq.gz>
//
// Build:
//
//	go mod init preproc_barcodes
//	go get github.com/biogo/hts@latest
//	go get github.com/klauspost/pgzip
//	go build -ldflags="-s -w" -o preproc_barcodes preproc_barcodes.go

import (
	"bytes"
	"regexp"

	"github.com/biogo/hts/sam"
)

var HaploRe = regexp.MustCompile(`(?:[ABCD][0-9]{2}){4}`)
var StlfrRe = regexp.MustCompile("([0-9]+_[0-9]+_[0-9]+)$")
var TellseqRe = regexp.MustCompile("([ATCGN]+)$")
var InvalidRe = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")

// ── SAM formatting ────────────────────────────────────────────────────────────

// pairedFlag returns the SAM FLAG integer for a paired unmapped read.
// R1: 0x1|0x4|0x8|0x40 = 77   R2: 0x1|0x4|0x8|0x80 = 141
func PairedFlag(isRead1 bool) int {
	if isRead1 {
		return 77
	}
	return 141
}

// workerState holds reusable per-goroutine scratch buffers.
type WorkerState struct {
	seqBuf  []byte       // assembled seq (pad + barcode + biological)
	qualBuf []byte       // assembled qual in ASCII Phred+33
	outBuf  bytes.Buffer // SAM text accumulator for the whole batch
}

// writeSAMRecord appends a single SAM record as a text line to ws.outBuf.
// qual is expected in ASCII Phred+33 (as read from FASTQ).
// This avoids constructing a sam.Record object entirely for the write path.
func (ws *WorkerState) writeuSAMRecord(name string, flag int, seq, qual []byte) {
	ws.outBuf.WriteString(name)
	ws.outBuf.WriteByte('\t')
	// write flag as decimal
	ws.writeInt(flag)
	// RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN — all * or 0 for uSAM
	ws.outBuf.WriteString("\t*\t0\t255\t*\t*\t0\t0\t")
	ws.outBuf.Write(seq)
	ws.outBuf.WriteByte('\t')
	ws.outBuf.Write(qual)
	ws.outBuf.WriteByte('\n')
}

// writeInt writes a non-negative integer to outBuf without allocating.
func (ws *WorkerState) writeInt(n int) {
	if n == 0 {
		ws.outBuf.WriteByte('0')
		return
	}
	var tmp [10]byte
	i := len(tmp)
	for n > 0 {
		i--
		tmp[i] = byte('0' + n%10)
		n /= 10
	}
	ws.outBuf.Write(tmp[i:])
}

// ─── barcode information ───────────────────────────────────────────────────
// Idenftify the barcode type as one of haplotagging, stlfr, or tellseq
// and whether it is valid or not. Returns barcode type (string) and valid (bool)
func DetermineBarcode(bc string) (string, bool) {
	var barcode string

	invalid := InvalidRe.MatchString(bc)

	if HaploRe.MatchString(bc) {
		barcode = "haplotagging"
	} else if StlfrRe.MatchString(bc) {
		barcode = "stlfr"
	} else if TellseqRe.MatchString(bc) {
		barcode = "tellseq"
	}
	return barcode, !invalid
}

// ─── SAM tag helpers ──────────────────────────────────────────────────────────

func GetStringTag(r *sam.Record, tag string) (string, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			if s, ok := aux.Value().(string); ok {
				return s, true
			}
		}
	}
	return "", false
}

func GetIntTag(r *sam.Record, tag string) (int, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			if s, ok := aux.Value().(int); ok {
				return s, true
			}
		}
	}
	return 0, false
}
