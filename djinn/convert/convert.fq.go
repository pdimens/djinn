package convert

import (
	"bufio"
	"djinn/barcodes"
	"djinn/fastq"
	"fmt"
	"io"
	"os"
	"strconv"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func ConvertFq(fqs []string, convTo, prefix, bcmap string, threads int) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var bcs barcodes.Generator
	var converter func(*fastq.CoreFq, *bufio.Writer)

	// --- assess the data type and get the parser -----------------------------
	coreParser, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		if bcmap == "" {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s. If this is 10X data, a barcode file must be provided to identify inline barcodes.", fqs[0])
		} else {
			// SWAP TO THE 10X CONVERTER and return
			//fmt.Fprintln(os.Stderr, "Format detected: 10X")
			return ConvertFrom10X(fqs, convTo, prefix, bcmap, threads)
		}
	}

	// --- create the barcode conversion generator -------------------------------------
	switch convTo {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
		converter = (*fastq.CoreFq).ToHaplotagging
	case "stlfr":
		bcs = barcodes.NewStlfr()
		converter = (*fastq.CoreFq).ToStlfr
	case "tellseq":
		bcs = barcodes.NewTellseq()
		converter = (*fastq.CoreFq).ToTellseq
	case "10x":
		bcs = barcodes.NewTenX()
		converter = (*fastq.CoreFq).ToTenX
	default:
		return fmt.Errorf("unknown barcode type %q", convTo)
	}
	defer bcs.Close()

	seen := make(map[string][]byte, 4_000_000)
	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer

	// ── open reader ───────────────────────────────────────────────────────────
	for idx, i := range fqs { // iterate over files
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

		// iterate through records
		for {
			_rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			rec, ok := coreParser(_rec)
			// if no barcode was found, add invalid BX, put to write and continue
			// VX is already invalid when no barcode found, doesn't need to be explicitly declared
			if !ok {
				rec.Barcode = bcs.GetInvalid()
				converter(&rec, outfq.Writer)
				continue
			}

			if convertedBC, isPresent := seen[string(rec.Barcode)]; isPresent {
				// barcode previously seen, pull the converted barcode
				rec.Barcode = convertedBC
			} else {
				// barcode not yet observed, generate new barcode, add it to map
				n, ok := bcs.NextInto(bcBuf)
				if !ok {
					return fmt.Errorf("too many unique barcodes for the conversion technology requested — unable to generate more barcodes.")
				}
				// independent copy; bcBuf gets reused/overwritten next iteration
				newBC := append([]byte(nil), bcBuf[:n]...)
				seen[string(rec.Barcode)] = newBC
				rec.Barcode = newBC
			}
			// write to output buffer
			converter(&rec, outfq.Writer)
		}
		fqReader.Close()
		outfq.Close()
	}

	// ------ write barcodes to map ------------------------
	f, err := os.Create(prefix + ".bc.map")
	if err != nil {
		return err
	}
	defer f.Close()

	writer := bufio.NewWriter(f)
	defer writer.Flush()

	for key, val := range seen {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.Write(val)
		writer.WriteByte('\n')
	}

	return nil
}

func ConvertFrom10X(fqs []string, convTo, prefix, bcmap string, threads int) error {
	return nil
}
