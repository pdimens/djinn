package barcodes

import (
	"bufio"
	"bytes"
	_ "embed"
	"fmt"

	"github.com/klauspost/compress/zstd"
)

// The real 10X barcode whitelist, one barcode per line, zstd-compressed and
// embedded directly into the binary. See data/README.md for how this file
// is produced/updated.
//
//go:embed data/10x_whitelist.txt.zst
var tenXWhitelistZst []byte

const tenXListLen = 16

var invalidTenXList = func() (b [tenXListLen]byte) {
	for i := range b {
		b[i] = 'N'
	}
	return
}()

// TenXList is a Generator that streams the real, curated 10X barcode
// whitelist line by line (via the embedded, compressed asset), instead of
// enumerating every possible 16-mer the way Nucleotides/NewTenX does. Use
// this when you need actual 10X barcodes rather than an arbitrary ACGT
// sequence of the right length.
type TenXList struct {
	dec     *zstd.Decoder
	scanner *bufio.Scanner
}

// NewTenXList opens a fresh stream over the embedded 10X barcode whitelist.
// Call Close when done to release the decoder.
func NewTenXList() (*TenXList, error) {
	dec, err := zstd.NewReader(bytes.NewReader(tenXWhitelistZst))
	if err != nil {
		return nil, fmt.Errorf("barcodes: opening embedded 10X whitelist: %w", err)
	}
	return &TenXList{dec: dec, scanner: bufio.NewScanner(dec)}, nil
}

// NextInto copies the next barcode in the whitelist into dst and reports
// whether one was available. As with every other Generator, dst must be
// used/consumed before the next call: bufio.Scanner's Bytes() (what this
// copies from) is only valid until the next Scan(), so there is nothing
// left to alias past this call either way.
func (t *TenXList) NextInto(dst []byte) (int, bool) {
	if !t.scanner.Scan() {
		return 0, false
	}
	return copy(dst, t.scanner.Bytes()), true
}

func (t *TenXList) GetInvalid() []byte { return invalidTenXList[:] }
func (t *TenXList) MaxLen() int        { return tenXListLen }
func (t *TenXList) Close()             { t.dec.Close() }
