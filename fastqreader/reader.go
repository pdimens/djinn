// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package fastqreader

import (
	"bufio"
	"log"
	"strings"
)

/*
 * This structure represents a single read from a fastq file
 */
type FastQRecord struct {
	Read1     []byte
	ReadQual1 []byte
	Read2     []byte
	ReadQual2 []byte
	Barcode   []byte
	Valid     bool
	Header    string
}

/*
 * This struture reprensets a "fastQ" reader that can pull single records
 * as well as sets of records (on the same barcode) from a fastq file
 */
type FastQReader struct {
	Line          int
	LastBarcode   []byte
	DefferedError error
	Pending       *FastQRecord
	R1Source      *ZipReader
	R1Buffer      *bufio.Reader
	R2Source      *ZipReader
	R2Buffer      *bufio.Reader
}

/* Open a new fastQ file */
func OpenFastQ(R1 string, R2 string) (*FastQReader, error) {

	var res = new(FastQReader)
	var err error

	res.R1Source, err = FastZipReader(R1)
	if err != nil {
		return nil, err
	}
	res.R2Source, err = FastZipReader(R2)
	if err != nil {
		return nil, err
	}

	res.R1Buffer = bufio.NewReader(res.R1Source)
	res.R2Buffer = bufio.NewReader(res.R2Source)
	res.Line = 0
	return res, nil
}

/* Parse a read header and find the barcode. Return the sanitized header, barcode, and 1/0 whether it's valid or not */
func ParseHeader(seq_id, lrType string) (string, []byte, bool) {
	var _barcode []byte
	var _valid bool
	var bxMatches []string
	var vxMatches []string

	// stringify the input
	// get the first part before the whitespaces
	_id := strings.Fields(seq_id)[0]
	_header := _id[:len(_id)-2]

	// find the barcode based on the type
	if lrType == "standard" {
		bxMatches = standardRe.FindStringSubmatch(seq_id)
	} else if lrType == "haplotagging" {
		bxMatches = haplotaggingRe.FindStringSubmatch(seq_id)
	} else if lrType == "stlfr" {
		bxMatches = stlfrRe.FindStringSubmatch(_id)
	} else {
		bxMatches = tellseqRe.FindStringSubmatch(_id)
	}

	if len(bxMatches) <= 1 {
		return "", []byte(""), false
	}
	barcode := bxMatches[1]

	// barcode exists so we need to find if it's valid or not
	if lrType == "standard" {
		_valid = strings.Contains(seq_id, "VX:i:1")
	} else {
		vxMatches = invalidRe.FindStringSubmatch(barcode)
		if len(vxMatches) > 1 {
			_valid = false
		} else {
			_valid = true
		}
	}

	_barcode = []byte(barcode)

	return _header, _barcode, _valid
}

// Read a single record from a fastQ file
func (fqr *FastQReader) ReadOneLine(result *FastQRecord, lrType string) error {
	/* Search for the next start-of-record.*/
	for {
		fqr.Line++
		R1_line, err := fqr.R1Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		R2_line, err := fqr.R2Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		if R1_line[0] == byte('@') {
			/* Found it! */
			result.Header, result.Barcode, result.Valid = ParseHeader(string(R1_line[1:]), lrType)
			break
		} else {
			log.Printf("Bad line in R1: %v at %v", string(R1_line), fqr.Line)
			log.Printf("Bad line in R2: %v at %v", string(R2_line), fqr.Line)
		}
	}

	/* Load the 4 lines for this record */
	var fastq_lines [6][]byte
	for i := range 4 {
		// skip the line with the + sign
		if i == 2 {
			continue
		}
		var err error
		var line []byte
		line, err = fqr.R1Buffer.ReadBytes(byte('\n'))
		fastq_lines[i] = line[0 : len(line)-1]
		if err != nil {
			return err
		}
		line, err = fqr.R2Buffer.ReadBytes(byte('\n'))
		fastq_lines[i+3] = line[0 : len(line)-1]
		if err != nil {
			return err
		}
	}

	/* Assign them to the right fields in the FastQRecord struct */
	result.Read1 = fastq_lines[1]
	result.ReadQual1 = fastq_lines[2]
	result.Read2 = fastq_lines[4]
	result.ReadQual2 = fastq_lines[5]
	// MAYBE A THING FOR COMMENTS?

	return nil
}

// Read a single record from a fastQ file, intended for inferring the linked-read barcode type
func (fqr *FastQReader) InferReadType() (string, error) {
	var barcode_type string

	/* Search for the next start-of-record.*/
	for {
		fqr.Line++
		R1_line, err := fqr.R1Buffer.ReadString(byte('\n'))
		if err != nil {
			return "unknown", err
		}
		R2_line, err := fqr.R2Buffer.ReadString(byte('\n'))
		if err != nil {
			return "unknown", err
		}
		if R1_line[0] == byte('@') {
			/* Found it! */
			barcode_type = findBC(R1_line)
			break
		} else {
			log.Printf("Bad line in R1: %v at %v", string(R1_line), fqr.Line)
			log.Printf("Bad line in R2: %v at %v", string(R2_line), fqr.Line)
		}
	}
	return barcode_type, nil
}

func findBC(seq_id string) string {
	var bxMatches []string
	var _id string

	// stringify the input
	// get the first part before the whitespaces
	_id = strings.Fields(seq_id)[0]

	// find the barcode based on the type
	// standard
	bxMatches = standardRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return "standard"
	}
	bxMatches = haplotaggingRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return "haplotagging"
	}
	bxMatches = stlfrRe.FindStringSubmatch(_id)
	if len(bxMatches) > 1 {
		return "stlfr"
	}
	bxMatches = tellseqRe.FindStringSubmatch(_id)
	if len(bxMatches) > 1 {
		return "tellseq"
	}
	return "unknown"

}
