package convert

import (
	"bufio"
	"djinn/barcodes"
	"djinn/fastq"
	"fmt"
	"io"
	"os"
	"slices"
	"strconv"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func ConvertFq(fqs []string, convTo, prefix, bcmap string, threads int) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var err error
	var bcs barcodes.Generator
	var converter func(*fastq.CoreFq, *bufio.Writer)
	var fmtType string

	// --- assess the data type and get the parser -----------------------------
	coreParser, fmtType, err := fastq.CoreFqParser(fqs[0])
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

	if fmtType == convTo {
		return fmt.Errorf("trying to convert to and from identical format (%s)", convTo)
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
		bcs, err = barcodes.NewTenXList()
		if err != nil {
			return err
		}
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
	// ── init inventory and generator ──────────────────────────────────────────
	var err error
	var bcs barcodes.Generator
	var converter func(*fastq.CoreFq, *bufio.Writer)

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
		return fmt.Errorf("trying to convert to and from identical format (10X)")
	}
	defer bcs.Close()

	seen := make(map[string][]byte, 4_000_000)
	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer
	bclist, err := barcodes.LoadTenXSet()
	if err != nil {
		return err
	}

	// ---- FQ readers -------------------------
	R1Reader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
	}
	defer R1Reader.Close()

	R2Reader, err := fastx.NewReader(seq.DNA, fqs[1], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[1], err)
	}
	defer R2Reader.Close()

	// ---- FQ writers -------------------------
	R1Writer, err := xopen.Wopen(prefix + ".R1.fq.gz")
	if err != nil {
		return err
	}
	defer R1Writer.Close()

	R2Writer, err := xopen.Wopen(prefix + ".R1.fq.gz")
	if err != nil {
		return err
	}
	defer R2Writer.Close()

	var r1core fastq.CoreFq
	var r2core fastq.CoreFq
	var ok bool

	// iterate through records
Loop:
	for {
		// process R1, which doesn't need an R2
		_r1, r1err := R1Reader.Read()
		switch r1err {
		case io.EOF:

		case nil:
			r1core, ok = fastq.TenXR12Corefq(_r1, bclist)
			// if no barcode was found, add invalid BX, put to write and continue
			// VX is already invalid when no barcode found, doesn't need to be explicitly declared
			if !ok {
				r1core.Barcode = bcs.GetInvalid()
			}
			if convertedBC, isPresent := seen[string(r1core.Barcode)]; isPresent {
				// barcode previously seen, pull the converted barcode
				r1core.Barcode = convertedBC
			} else {
				// barcode not yet observed, generate new barcode, add it to map
				n, ok := bcs.NextInto(bcBuf)
				if !ok {
					return fmt.Errorf("too many unique barcodes for the conversion technology requested — unable to generate more barcodes.")
				}
				// independent copy; bcBuf gets reused/overwritten next iteration
				newBC := append([]byte(nil), bcBuf[:n]...)
				seen[string(r1core.Barcode)] = newBC
				r1core.Barcode = newBC
			}
			converter(&r1core, R1Writer.Writer)
		default:
			return err
		}

		// process R2
		_r2, r2err := R2Reader.Read()
		switch r2err {
		case io.EOF:
			// break when both are exhausted
			if r1err == io.EOF {
				break Loop
			}
		case nil:
			r2core = fastq.TenXR22Corefq(_r2)
			// check if R1 exists and match the seq ID
			if r1err == nil {
				if ok := slices.Equal(r1core.ID, r2core.ID); !ok {
					return fmt.Errorf("Misaligned R1 and R2 records. The 10X format only stores the barcode at the beginning of read 1, meaning Read 2 must have the same sequence ID as Read 1 to reliably assign a barcode to Read 2, but they do not match.\nRead 1:\n%s\nRead 2:\n%s", r1core.ID, r2core.ID)
				}
				r2core.Barcode = r1core.Barcode
				r2core.Valid = r1core.Valid
			}
			converter(&r2core, R2Writer.Writer)
		default:
			return err
		}
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
