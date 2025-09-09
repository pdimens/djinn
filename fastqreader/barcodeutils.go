package fastqreader

import (
	"log"
	"regexp"
)

var bxRe = regexp.MustCompile(`BX:Z:(\S+)\s`)
var vxRe = regexp.MustCompile(`VX:i:([01])\s`)

var haplotaggingRe = regexp.MustCompile(`BX:Z:(A\d{2}C\d{2}B\d{2}D\d{2})\s`)

var stlfrRe = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)\s`)
var stlfrInvalidRe = regexp.MustCompile(`^0_|_0_|_0$`)

var tellseqRe = regexp.MustCompile(`:([ATCGN]+)\s`)

/*Return true if the fastq record is in Standard format*/
func isStandardized(seq_id string) bool {
	// regex match BX:Z:*
	bxMatches := bxRe.FindStringSubmatch(seq_id)
	if len(bxMatches) <= 1 {
		return false
	}
	// regex match VX:i:[01]
	vxMatches := vxRe.FindStringSubmatch(seq_id)
	if len(vxMatches) > 1 {
		return true
	} else {
		return false
	}
}

/*Return true if the fastq record is in haplotagging format*/
func isHaplotagging(seq_id string) bool {
	// regex match BX:Z:AxxCxxBxxDxx
	bxMatches := haplotaggingRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return true
	}
	return false
}

/*Return true if the fastq record is in stlfr format*/
func isStlfr(seq_id string) bool {
	// regex match #1_2_3
	bxMatches := stlfrRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return true
	}
	return false
}

/*Return true if the fastq record is in stlfr format*/
func isTellseq(seq_id string) bool {
	// regex match :ATCGN
	bxMatches := tellseqRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return true
	}
	return false
}

/*
Parse through the first 200 records of a paired-end fastq and figure out what kind of
format it's in. Returns early if a format is detected.
*/
func FindFastqFormat(r1, r2 string) (string, error) {
	var record FastQRecord
	var err error
	fqr, err := OpenFastQ(r1, r2)

	if err != nil {
		return "", err
	}

	for range 200 {
		err := fqr.ReadOneLine(&record)
		if err != nil {
			return "", err
		}
		if isStandardized(record.Header) {
			return "standard", nil
		} else if isHaplotagging(record.Header) {
			return "haplotagging", nil
		} else if isStlfr(record.Header) {
			return "stlfr", nil
		} else if isTellseq(record.Header) {
			return "tellseq", nil
		}
	}
	log.Fatal("Unable to identify the input FASTQ files as one of haplotagging, stLFR, or TELLseq. Unable to proceed.")
	return "unknown", nil
}
