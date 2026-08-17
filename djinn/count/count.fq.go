package count

import (
	"bufio"
	"djinn/fastq"
	"fmt"
	"io"
	"maps"
	"os"
	"strconv"

	"github.com/shenwei356/bio/seqio/fastx"
)

func CountFQ(fqs []string, invalid bool) error {
	// determine what kind of linked-read tech it is
	processBC, err := fastq.CheckFastqFormat(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}

	// ── open barcode writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	// ---- main loop -------------------------
	set := make(map[string]int, 7_000_000) // barcode container
	var bc string
	var valid bool

	// R1
	fqReader, err := fastx.NewDefaultReader(fqs[0])
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
	}

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

		set[bc]++
	}
	fqReader.Close()

	// R2
	if len(fqs) > 1 {
		setR2 := make(map[string]int, 500_000) // R2 barcode container

		fqReader, err = fastx.NewDefaultReader(fqs[1])
		if err != nil {
			return fmt.Errorf("opening %s: %w", fqs[1], err)
		}

		for {
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

			// if the barcode already appears in the R1 set, skip
			if _, ok := set[bc]; ok {
				continue
			}
			setR2[bc]++
		}
		maps.Copy(set, setR2)
	}

	for key, val := range set {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.WriteString(strconv.Itoa(int(val)))
		writer.WriteByte('\n')
		fmt.Printf("%s\t%d\n", key, val)
	}
	return nil
}
