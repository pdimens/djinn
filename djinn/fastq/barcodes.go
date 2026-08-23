package fastq

import (
	"bytes"
	"fmt"

	regexp "github.com/coregx/coregex"
	"github.com/shenwei356/bio/seqio/fastx"
)

var MissingBarcode []byte

var StlfTell = regexp.MustCompile(`(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`)
var BxBarcode = regexp.MustCompile(`(?:BX\:Z\:(\S+))`)
var Invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")
var Tellseq = regexp.MustCompile(`:([ATCGN]+)(\s|$)`)

var Stlfr = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)(\s|$)`)
var StdBx = regexp.MustCompile(`BX:Z:(\S+)(?:\s|$)`)
var StdVx = regexp.MustCompile(`VX:i:([01])(?:\s|$)`)
var IlluminaOld = regexp.MustCompile(`/[12](?:\s|$)`)
var IlluminaNew = regexp.MustCompile(`[12]:[YN]:\d+:[A-Za-z0-9]+(?:\s|$)`)

func FindBarcode(rec *fastx.Record) string {
	var bxVal []byte
	bxMatch := BxBarcode.FindSubmatch(rec.Desc)
	if bxMatch == nil {
		matches := StlfTell.FindSubmatch(rec.ID)
		if matches != nil {
			switch {
			case len(matches) > 1 && matches[1] != nil:
				// matches[1] is the tellseq barcode e.g. "ATCGN"
				bxVal = matches[1]
			case len(matches) > 2 && matches[2] != nil:
				// matches[2] is the stlfr barcode e.g. "1_2_3"
				bxVal = matches[2]
			}
		}
	} else {
		bxVal = bxMatch[1]
	}
	return string(bxVal)
}

// Detect the linked-read technology type from the first 100 records of the FASTQ file.
// Returns the function to be used to detect barcodes and process reads in all records within the main loop.
func CheckFastqFormat(fq string) (func(rec *fastx.Record) (string, bool), error) {
	var rec *fastx.Record
	var h, t, s int
	var totalReads int

	fqReader, err := fastx.NewDefaultReader(fq)
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", fq, err)
	}
	defer fqReader.Close()

	for i := range 100 {
		rec, err = fqReader.Read()
		if err != nil {
			return nil, fmt.Errorf("reading %v, record %v: %w", fq, i, err)
		}
		// is there a BX tag in the comments/description?
		if StdBx.Match(rec.Desc) {
			h += 1
		}
		// if not, look for tellseq
		if Tellseq.Match(rec.ID) {
			t += 1
		}
		// if not, look for stlfr
		if Stlfr.Match(rec.ID) {
			s += 1
		}
		totalReads += 1
	}
	// if more than one style found, return an error
	// otherwise, set global missing barcode for that chemistry and parsing/standardizing function
	if (h + s + t) > totalReads {
		return nil, fmt.Errorf("more than one linked-read technology format identified. Input data must use a single format. Reads types identified: Haplotagging - %d | stLFR - %d | TELLseq - %d.", h, s, t)
	} else if h > 0 {
		MissingBarcode = []byte("VX:i:0\tBX:Z:A00C00B00D00")
		return FqHaplotagBX, nil
	} else if s > 0 {
		MissingBarcode = []byte("VX:i:0\tBX:Z:0_0_0")
		return FqStlfr, nil
	} else if t > 0 {
		MissingBarcode = []byte("VX:i:0\tBX:Z:NNNNNNNNNNNNNNNNNN")
		return FqTellseq, nil
	} else {
		return nil, fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fq)
	}
}

func FqHaplotagBX(rec *fastx.Record) (string, bool) {
	var valid bool
	var bc []byte

	bloc := StdBx.FindSubmatchIndex(rec.Desc)
	if bloc == nil {
		return "", false
	}

	bc = rec.Desc[bloc[2]:bloc[3]]
	vloc := StdVx.FindSubmatchIndex(rec.Desc)
	if vloc == nil {
		// no VX tag, check if barcode valid
		if !Invalid.Match(bc) {
			valid = true
		}
	} else if bytes.Equal(rec.Desc[vloc[2]:vloc[3]], []byte{'1'}) {
		// VX tag present and is 1
		valid = true
	}

	return string(bc), valid
}

// Find tellseq barcode inline in rec.ID, checking it for validity.
// Converts rec in-place into standard format. Returns early if no barcode was found.
// Returns a bool of whether the barcode was valid (true) or not (false) as a sentinel
// value for how to post-process the read.
func FqTellseq(rec *fastx.Record) (string, bool) {
	var valid bool
	var bc []byte

	// find stlfr barcode in the record ID
	bloc := Tellseq.FindSubmatchIndex(rec.ID)
	if bloc == nil {
		return "", false
	}

	bc = rec.ID[bloc[2]:bloc[3]]
	if !Invalid.Match(bc) {
		valid = true
	}

	return string(bc), valid
}

// Find tellseq barcode inline in rec.ID, checking it for validity.
// Converts rec in-place into standard format. Returns early if no barcode was found.
// Returns a bool of whether the barcode was valid (true) or not (false) as a sentinel
// value for how to post-process the read.
func FqStlfr(rec *fastx.Record) (string, bool) {
	var valid bool
	var bc []byte

	// find stlfr barcode
	bloc := Stlfr.FindSubmatchIndex(rec.ID)
	if bloc == nil {
		return "", false
	}

	bc = rec.ID[bloc[2]:bloc[3]]
	if !Invalid.Match(bc) {
		valid = true
	}

	return string(bc), valid
}
