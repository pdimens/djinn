package extract

import (
	"bufio"
	"djinn/fastq"
	"fmt"
	"io"
	"os"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func ExtractFQ(fqs []string, invalid bool) error {
	// determine what kind of linked-read tech it is
	processBC, err := fastq.CheckFastqFormat(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	seq.ValidateSeq = false
	// ── open barcode writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	// ---- main loop -------------------------
	set := make(map[string]struct{}, 7_000_000) // barcode container

	for _, fq := range fqs {
		if err := extractFromFastq(fq, invalid, processBC, set, writer); err != nil {
			return err
		}
	}
	return nil
}

// extractFromFastq reads a single FASTQ file and writes newly-seen barcodes
// to writer. It closes its reader before returning, rather than deferring
// the close to the caller's lifetime (which, across many input files, would
// keep every prior file handle open until ExtractFQ itself returned).
func extractFromFastq(fq string, invalid bool, processBC func(*fastx.Record) (string, bool), set map[string]struct{}, writer *bufio.Writer) error {
	fqReader, err := fastx.NewReader(seq.DNA, fq, "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fq, err)
	}
	defer fqReader.Close()

	var bc string
	var valid bool
	for {
		// iterate through records
		rec, err := fqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		bc, valid = processBC(rec)
		if bc == "" {
			continue
		} // no barcode

		if !valid && !invalid {
			continue
		} // is invalid and was not asked to report invalid

		if _, ok := set[bc]; ok {
			continue
		} else {
			set[bc] = struct{}{}
			writer.WriteString(bc)
			writer.WriteByte('\n')
		}
	}
	return nil
}
