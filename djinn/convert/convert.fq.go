package convert

import (
	"bufio"
	"bytes"
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

func ConvertFq(fqs []string, convTo, prefix string, threads int) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var bcs barcodes.Generator
	var converter func(*fastq.CoreFq, *bytes.Buffer)
	//var fmtConverter func(rec *fastq.CoreFq, writer *bytes.Buffer)

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

	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer

	// ── open reader ───────────────────────────────────────────────────────────
	for idx, i := range fqs { // iterate over files
		// determine what kind of linked-read tech it is
		lrFormat, err := fastq.CheckLRFormat(fqs[0])
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

		for { // iterate through records
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

		}
		fqReader.Close()
		outfq.Close()
	}

	seen := make(map[string][]byte, 4_000_000)

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
