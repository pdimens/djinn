package ncbi

import (
	"bufio"
	"djinn/barcodes"
	"djinn/fastq"
	"fmt"
	"io"
	"os"
	"slices"

	"github.com/biogo/hts/bgzf"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func NcbiFq(fqs []string, bcmap string, threads int) error {
	var err error

	// --- assess the data type and get the parser -----------------------------
	coreParser, _, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		if bcmap == "" {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s. If this is 10X data, a barcode file must be provided to identify inline barcodes.", fqs[0])
		} else {
			// SWAP TO THE 10X CONVERTER and return
			return NCBIfq10x(fqs, bcmap, threads)
		}
	}

	// ── open reader ───────────────────────────────────────────────────────────
	for _, i := range fqs { // iterate over files
		// ---- FQ reader -------------------------
		fqReader, err := fastx.NewReader(seq.DNA, i, "")
		if err != nil {
			return fmt.Errorf("opening %s: %w", i, err)
		}

		// ---- BAM writer -------------------------
		bgzfW := bgzf.NewWriter(os.Stdout, threads) // concurrent blocks
		bw := bufio.NewWriter(bgzfW)

		// ---- iterate through records -------------------------
		for {
			_rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			rec, _ := coreParser(_rec)
			rec.WriteBam(bw)
		}
		fqReader.Close()
	}
	if err := os.Stdout.Close(); err != nil {
		return err
	}
	return nil
}

func NCBIfq10x(fqs []string, bcmap string, threads int) error {
	var err error

	// ---- read in barcodes -------------------------------
	bclist, err := barcodes.LoadTenXFromFile(bcmap)
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

	// ---- BAM writer -------------------------
	bgzfW := bgzf.NewWriter(os.Stdout, threads)
	bw := bufio.NewWriter(bgzfW)
	bw.WriteString("@HD	VN:1.6	SO:unsorted\n")

	var r1core fastq.CoreFq
	var r2core fastq.CoreFq

Loop:
	for {
		// process R1, which doesn't need an R2
		_r1, r1err := R1Reader.Read()
		switch r1err {
		case io.EOF:

		case nil:
			r1core, _ = fastq.TenXR12Corefq(_r1, bclist)
			r1core.WriteBam(bw)
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
			r2core.WriteBam(bw)
		default:
			return err
		}
	}

	if err := os.Stdout.Close(); err != nil {
		return err
	}
	return nil
}
