package barcodes

import (
	"bufio"
	"bytes"
	"compress/gzip"
	_ "embed"
	"fmt"
)

// The real 10X Genomics barcode whitelist ("4M-with-alts-february-2016"):
// one 16bp barcode per line, gzip-compressed, embedded directly into the
// binary. 4,792,320 barcodes, all exactly 16 bases.
//
//go:embed 10x.gz
var tenXGz []byte

const tenXLen = 16

var invalidTenX = func() (b [tenXLen]byte) {
	for i := range b {
		b[i] = 'N'
	}
	return
}()

// TenXList is a Generator that streams the real, curated 10X barcode
// whitelist line by line (from the embedded asset), instead of
// enumerating every possible 16-mer the way Nucleotides/NewTenX does. Use
// this when assigning *new* 10X barcodes during a conversion -- e.g. the
// "10x" case in convert.fq.go's output-format switch.
type TenXList struct {
	gz      *gzip.Reader
	scanner *bufio.Scanner
}

// NewTenXList opens a fresh stream over the embedded 10X barcode whitelist.
// Call Close when done to release the decoder.
func NewTenXList() (*TenXList, error) {
	gz, err := gzip.NewReader(bytes.NewReader(tenXGz))
	if err != nil {
		return nil, fmt.Errorf("barcodes: opening embedded 10X whitelist: %w", err)
	}
	return &TenXList{gz: gz, scanner: bufio.NewScanner(gz)}, nil
}

// NextInto copies the next barcode in the whitelist into dst and reports
// whether one was available. As with every other Generator, dst must be
// used/consumed before the next call: bufio.Scanner's Bytes() (what this
// copies from) is only valid until the next Scan() anyway.
func (t *TenXList) NextInto(dst []byte) (int, bool) {
	if !t.scanner.Scan() {
		return 0, false
	}
	return copy(dst, t.scanner.Bytes()), true
}

func (t *TenXList) GetInvalid() []byte { return invalidTenX[:] }
func (t *TenXList) MaxLen() int        { return tenXLen }
func (t *TenXList) Close()             { t.gz.Close() }

// LoadTenXSet reads the full embedded 10X barcode whitelist into a
// membership set, for recognizing/validating an inline 10X barcode parsed
// out of a read's sequence (see fastq.TenX2Corefq's bc_map parameter) --
// a different access pattern than TenXList's sequential NextInto: this one
// needs the whole whitelist resident as a set for O(1) lookups, rather
// than handing out barcodes one at a time.
func LoadTenXSet() (map[string]struct{}, error) {
	gz, err := gzip.NewReader(bytes.NewReader(tenXGz))
	if err != nil {
		return nil, fmt.Errorf("barcodes: opening embedded 10X whitelist: %w", err)
	}
	defer gz.Close()

	set := make(map[string]struct{}, 5_000_000)
	scanner := bufio.NewScanner(gz)
	for scanner.Scan() {
		// map keys always copy on insert (Go strings are immutable and a
		// []byte->string conversion here copies), so aliasing
		// scanner.Bytes() past this call is a non-issue.
		set[scanner.Text()] = struct{}{}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("barcodes: reading embedded 10X whitelist: %w", err)
	}
	return set, nil
}
