package invalid

import (
	"djinn/fastq"
	"fmt"
	"io"
	"strconv"

	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func FilterInvalidFQ(fqs []string, prefix, invalidprefix string) error {
	keepInvalid := invalidprefix != ""

	for idx, i := range fqs { // iterate over files
		// determine what kind of linked-read tech it is
		processBC, err := fastq.CheckFastqFormat(fqs[0])
		if err != nil {
			return fmt.Errorf("%w", err)
		}
		// ---- FQ reader -------------------------
		fqReader, err := fastx.NewDefaultReader(i)
		if err != nil {
			return fmt.Errorf("opening %s: %w", i, err)
		}

		// ---- FQ writer -------------------------
		outfq, err := xopen.Wopen(prefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz")
		if err != nil {
			return err
		}

		var outfqInv *xopen.Writer
		if keepInvalid {
			outfqInv, err = xopen.Wopen(invalidprefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz")
			if err != nil {
				return err
			}
		}
		var bc string
		var valid bool

		for { // iterate through records
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

			bc, valid = processBC(rec)
			if bc == "" || !valid {
				if keepInvalid {
					rec.FormatToWriter(outfqInv, 0)
				}
				continue
			}

			rec.FormatToWriter(outfq, 0)
		}
		fqReader.Close()
		outfq.Close()
		if keepInvalid {
			outfqInv.Close()
		}
	}

	return nil
}
