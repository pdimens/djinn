package ncbi

import (
	"bufio"
	"djinn/fastq"
	"fmt"
	"io"
	"os"

	"github.com/biogo/hts/bgzf"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func NcbiFq(infiles []string, threads int, asSam bool) error {
	readThread := 1
	writeThread := 1
	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}

	// ── assess read and barcode style ────────────────────────────────
	processBC, err := fastq.CheckFastqFormat(infiles[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}

	// ── FASTQ readers ──────────────────────────────────────────────
	R1Reader, err := fastx.NewReader(seq.DNA, infiles[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", infiles[0], err)
	}
	defer R1Reader.Close()

	R2Reader, err := fastx.NewReader(seq.DNA, infiles[1], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", infiles[1], err)
	}
	defer R2Reader.Close()

	// ── open XAM writer ───────────────────────────────────────────────────────────

	bgzfW := bgzf.NewWriter(os.Stdout, writeThread) // concurrent blocks
	bw := bufio.NewWriter(bgzfW)

	// ── loop through fastqs  --------──────────────────────────────────────────────
	var bc string
	var valid bool

	for {
		// R1 record
		rec, err := R1Reader.Read()
		if err != nil && err != io.EOF { // read error and not EOF
			return err
		}
		if err != io.EOF { // process record if not EOF
			bc, valid = processBC(rec)
			// create sam record and write into channel
		}

		// R2 record
		rec2, err2 := R2Reader.Read()
		if err2 != nil && err2 != io.EOF { // read error and not EOF
			return err
		} else if err == io.EOF && err2 == io.EOF { // both files are EOF, we're done
			break
		}
		if err2 != io.EOF { // process record if not EOF
			bc, valid = processBC(rec2)
			// create sam record and write into channel
		}
	}

	return nil
}
