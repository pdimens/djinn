package fastq

import (
	"bytes"
	"fmt"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// Detect the linked-read technology type from the first 100 records of the FASTQ file.
// Returns the function to be used to detect barcodes and process reads in all records within the main loop.
func CheckFastqFormat(fq string) (func(rec *fastx.Record) (string, bool), error) {
	var rec *fastx.Record
	var h, t, s int
	var totalReads int
	seq.ValidateSeq = false

	fqReader, err := fastx.NewReader(seq.DNA, fq, "")
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
