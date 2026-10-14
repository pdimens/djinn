package xam

import (
	"bufio"
	"djinn/barcodes"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// isHaplotaggingBarcode reports whether bc has the fixed haplotagging
// shape "A##C##B##D##" -- 12 bytes: a letter followed by two digits,
// repeated for A, C, B, D in that order.
func isHaplotaggingBarcode(bc string) bool {
	if len(bc) != 12 {
		return false
	}
	letters := [4]byte{'A', 'C', 'B', 'D'}
	for i, want := range letters {
		pos := i * 3
		if bc[pos] != want || !isDigit(bc[pos+1]) || !isDigit(bc[pos+2]) {
			return false
		}
	}
	return true
}

// isNucleotideBarcode reports whether bc consists solely of DNA bases
// (including the ambiguity code N) -- the shape shared by both TELLseq
// and 10X barcodes once standardized into a BX:Z tag. Distinguishing the
// two requires checking against the 10X whitelist; see
// DetectBarcodeFormat.
func isNucleotideBarcode(bc string) bool {
	if bc == "" {
		return false
	}
	for i := 0; i < len(bc); i++ {
		if !isBase(bc[i]) {
			return false
		}
	}
	return true
}

// openAlignmentReader opens file (or stdin for "-") as a SAM or BAM
// stream, auto-detecting which from the first two bytes. It's meant for
// a short, throwaway read (see DetectBarcodeFormat): SAM/BAM streams
// can't be rewound, so a caller that goes on to actually process the
// file must open a fresh reader of its own afterward -- this one should
// just be closed when done sniffing.
//
// Closing isBAM's fh is the caller's responsibility when isBAM is true:
// *bam.Reader's own Close() only releases its internal decompression
// resources, not the underlying file handle (matching the asymmetry
// NewXamReaderChan's goroutine already works around the same way; SAM's
// wrapper Close() does close its file handle, via SamReader.closer).
func openAlignmentReader(file string) (r AlignmentReader, fh *os.File, isBAM bool, err error) {
	fh = os.Stdin
	if file != "-" {
		fh, err = os.Open(file)
		if err != nil {
			return nil, nil, false, err
		}
	}

	br := bufio.NewReaderSize(fh, IoBuf)
	magic, err := br.Peek(2)
	if err != nil {
		if fh != os.Stdin {
			fh.Close()
		}
		return nil, nil, false, err
	}

	isBAM = magic[0] == 0x1f && magic[1] == 0x8b
	if isBAM {
		br2, err := bam.NewReader(br, 1)
		if err != nil {
			if fh != os.Stdin {
				fh.Close()
			}
			return nil, nil, false, err
		}
		return br2, fh, true, nil
	}

	sr, err := sam.NewReader(br)
	if err != nil {
		if fh != os.Stdin {
			fh.Close()
		}
		return nil, nil, false, err
	}
	return &SamReader{Reader: sr, closer: fh}, fh, false, nil
}

// DetectBarcodeFormat scans up to the first 100 records of file (SAM or
// BAM, "-" for stdin) and reports which linked-read barcode format they
// use: "haplotagging", "stlfr", "tellseq", "10x", or "error" when the
// format can't be determined -- no recognizable barcode on any sampled
// record, the file can't be opened/read, or more than one format was
// identified (input is expected to be consistently one format,
// mirroring fastq.CoreFqParser's own rule).
//
// Unlike fastq.CoreFqParser, this can and does report "10x": a SAM/BAM
// barcode is always stored as a BX:Z tag regardless of its original
// source technology, so TELLseq and 10X -- both plain nucleotide
// strings once standardized -- are otherwise indistinguishable by shape
// alone. They're told apart by checking against the embedded 10X
// whitelist (barcodes.LoadTenXSet), loaded only if a nucleotide-shaped
// barcode is actually seen in the sample.
func DetectBarcodeFormat(file string) string {
	r, fh, isBAM, err := openAlignmentReader(file)
	if err != nil {
		return "error"
	}
	defer r.Close()
	if isBAM && fh != os.Stdin {
		defer fh.Close()
	}

	var haplotagging, stlfr int
	var nucBarcodes []string
	for i := 0; i < 100; i++ {
		rec, err := r.Read()
		if err != nil {
			break // EOF or a read error -- classify from whatever was sampled
		}
		bc, ok := FindBarcode(rec)
		if !ok || bc == "" {
			continue
		}
		switch {
		case isHaplotaggingBarcode(bc):
			haplotagging++
		case isStlfrBody(bc):
			stlfr++
		case isNucleotideBarcode(bc):
			nucBarcodes = append(nucBarcodes, bc)
		}
	}

	var tenx, tellseq int
	if len(nucBarcodes) > 0 {
		set, err := barcodes.LoadTenXSet()
		if err != nil {
			return "error"
		}
		for _, bc := range nucBarcodes {
			if _, ok := set[bc]; ok {
				tenx++
			} else {
				tellseq++
			}
		}
	}

	formatsFound := 0
	for _, n := range [4]int{haplotagging, stlfr, tenx, tellseq} {
		if n > 0 {
			formatsFound++
		}
	}

	switch {
	case formatsFound != 1:
		return "error"
	case haplotagging > 0:
		return "haplotagging"
	case stlfr > 0:
		return "stlfr"
	case tenx > 0:
		return "10x"
	default:
		return "tellseq"
	}
}
