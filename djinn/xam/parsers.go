package xam

import (
	"bufio"
	"fmt"
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
// and 10X barcodes once standardized into a BX:Z tag. The two are told
// apart purely by length (10X: 16, TELLseq: 18); see DetectBarcodeFormat.
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
// use: "haplotagging", "stlfr", "genetic-<LEN>", or "error" when the
// format can't be determined -- no recognizable barcode on any sampled
// record, the file can't be opened/read, or more than one format was
// identified (input is expected to be consistently one format,
// mirroring fastq.CoreFqParser's own rule).
//
// A SAM/BAM barcode is always stored as a BX:Z tag regardless of its
// original source technology, so TELLseq and 10X -- both plain
// nucleotide strings once standardized -- are otherwise
// indistinguishable by shape alone. Rather than resolving that via the
// (expensive to load) 10X whitelist, nucleotide-shaped barcodes are
// reported as "genetic-<LEN>" (e.g. "genetic-16" for 10X, "genetic-18"
// for TELLseq) and left for the caller to interpret.
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
	nucLens := make(map[int]int)
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
			nucLens[len(bc)]++
		}
	}

	formatsFound := 0
	if haplotagging > 0 {
		formatsFound++
	}
	if stlfr > 0 {
		formatsFound++
	}
	formatsFound += len(nucLens)

	if formatsFound != 1 {
		return "error"
	}
	switch {
	case haplotagging > 0:
		return "haplotagging"
	case stlfr > 0:
		return "stlfr"
	default:
		for length := range nucLens {
			return fmt.Sprintf("genetic-%d", length)
		}
		return "error" // unreachable: formatsFound==1 guarantees one of the above
	}
}
