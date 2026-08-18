package singletons

import (
	"djinn/fastq"
	"fmt"
	"io"
	"maps"

	"github.com/shenwei356/bio/seqio/fastx"
)

func getFqCount(infiles []string) (map[string]int16, error) {
	set := make(map[string]int16, 7_000_000) // barcode container
	var emptymap map[string]int16            // for returning empty thing on error
	var bc string
	var valid bool

	// determine what kind of linked-read tech it is
	processBC, err := fastq.CheckFastqFormat(infiles[0])
	if err != nil {
		return emptymap, err
	}

	// R1
	fqReader, err := fastx.NewDefaultReader(infiles[0])
	if err != nil {
		return emptymap, fmt.Errorf("opening %s: %w", infiles[0], err)
	}

	for {
		// iterate through records
		rec, err := fqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return emptymap, err
		}
		bc, valid = processBC(rec)
		if bc == "" || !valid {
			continue
		} // no barcode or is invalid
		set[bc]++
	}
	fqReader.Close()

	// R2
	if len(infiles) > 1 {
		setR2 := make(map[string]int16, 500_000) // R2 barcode container

		fqReader, err = fastx.NewDefaultReader(infiles[1])
		if err != nil {
			return emptymap, fmt.Errorf("opening %s: %w", infiles[1], err)
		}
		defer fqReader.Close()

		for {
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return emptymap, err
			}

			bc, valid = processBC(rec)
			if bc == "" || !valid {
				continue
			} // no barcode or is invalid
			if _, ok := set[bc]; ok {
				continue
			} // if the barcode already appears in the R1 set, skip
			setR2[bc]++
		}
		maps.Copy(set, setR2)
	}
	return set, nil
}

func FilterSingletons(infiles []string, singletons string) error {
	// guard against draining stdin when getting barcode counts
	bcCounts, err := getFqCount(infiles)
	if err != nil {
		return err
	}
	// ── open reader ───────────────────────────────────────────────────────────

	// ── open writer ───────────────────────────────────────────────────────────

	// ── loop record channel ──────────────────────────────────────────────

	// if bcCounts[bxVal] > 2 {
	// write record
	// }
	return nil
}
