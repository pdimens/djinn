package singletons

import (
	"bufio"
	"djinn/fastq"
	"fmt"
	"io"
	"maps"
	"os"
	"strconv"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func getFqCount(infiles []string) (map[string]int16, error) {
	set := make(map[string]int16, 7_000_000) // barcode container
	var emptymap map[string]int16            // for returning empty thing on error
	var bc string
	var valid bool
	seq.ValidateSeq = false

	// determine what kind of linked-read tech it is
	processBC, err := fastq.CheckFastqFormat(infiles[0])
	if err != nil {
		return emptymap, err
	}

	// R1
	fqReader, err := fastx.NewReader(seq.DNA, infiles[0], "")
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

func FilterSingletonsFQ(fqs []string, prefix, singletonprefix, barcodecount string) error {
	// guard against draining stdin when getting barcode counts
	bcCounts, err := getFqCount(fqs)
	if err != nil {
		return err
	}
	if barcodecount != "" {
		f, err := os.Create(barcodecount)
		if err != nil {
			return err
		}
		writer := bufio.NewWriter(f)
		for key, val := range bcCounts {
			writer.WriteString(key)
			writer.WriteByte('\t')
			writer.WriteString(strconv.Itoa(int(val)))
			writer.WriteByte('\n')
		}
		writer.Flush()
		f.Close()
	}
	keepSingle := singletonprefix != ""

	for idx, i := range fqs { // iterate over files
		// determine what kind of linked-read tech it is
		processBC, err := fastq.CheckFastqFormat(fqs[0])
		if err != nil {
			return fmt.Errorf("%w", err)
		}
		// ---- FQ reader -------------------------
		fqReader, err := fastx.NewReader(seq.DNA, i, "")
		if err != nil {
			return fmt.Errorf("opening %s: %w", i, err)
		}

		// ---- FQ writer -------------------------
		outfq, err := xopen.Wopen(prefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz")
		if err != nil {
			return err
		}

		var outfqSingle *xopen.Writer
		if keepSingle {
			outfqSingle, err = xopen.Wopen(singletonprefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz")
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
				continue
			}
			if bcCounts[bc] >= 2 {
				rec.FormatToWriter(outfq, 0)
			} else if keepSingle {
				rec.FormatToWriter(outfqSingle, 0)
			}
		}
		fqReader.Close()
		outfq.Close()
		if keepSingle {
			outfqSingle.Close()
		}
	}

	return nil
}
