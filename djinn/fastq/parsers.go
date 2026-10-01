package fastq

import (
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// Scan a fastq file and return the type of reading parser it needs, i.e., a haplotagging one,
// tellseq, stlfr, 10x
func CoreFqParser(fq string) (func(rec *fastx.Record) (core CoreFq, ok bool), error) {
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
			if errors.Is(err, io.EOF) {
				break
			}
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
	formatsFound := 0
	if h > 0 {
		formatsFound++
	}
	if s > 0 {
		formatsFound++
	}
	if t > 0 {
		formatsFound++
	}
	if formatsFound > 1 {
		return nil, fmt.Errorf("more than one linked-read technology format identified. Input data must use a single format. Reads types identified: Haplotagging - %d | stLFR - %d | TELLseq - %d.", h, s, t)
	}
	switch {
	case h > 0:
		fmt.Fprintln(os.Stderr, "Format detected: Haplotagging")
		return Haplotag2Corefq, nil
	case s > 0:
		fmt.Fprintln(os.Stderr, "Format detected: stLFR")
		return Stlfr2Corefq, nil
	case t > 0:
		fmt.Fprintln(os.Stderr, "Format detected: TELL-seq")
		return Tellseq2Corefq, nil
	default:
		return nil, nil
	}
}
