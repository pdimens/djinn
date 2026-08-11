package ncbi

import (
	"djinn/fastq"
	"djinn/xam"
	"fmt"

	"github.com/biogo/hts/sam"
)

func NCBI(infile, prefix string, threads int) error {
	r1Path := prefix + ".R1.fq.gz"
	r2Path := prefix + ".R2.fq.gz"

	readThread := 1
	writeThread := 1
	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}
	// ── open reader ───────────────────────────────────────────────────────────
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── FASTQ writers ──────────────────────────────────────────────
	threadsPerWriter := max(writeThread/2, 1)
	fw1, err := fastq.NewFastqWriter(r1Path, "/1", 4, threadsPerWriter)
	if err != nil {
		return fmt.Errorf("Error opening R1 output: %v\n", err)
	}
	defer fw1.Close()

	fw2, err := fastq.NewFastqWriter(r2Path, "/2", 4, threadsPerWriter)
	if err != nil {
		return fmt.Errorf("Error opening R2 output: %v\n", err)
	}
	defer fw2.Close()

	// ── loop record channel ──────────────────────────────────────────────
	var fw *fastq.FastqWriter

	for rec := range recChan {
		if rec.Flags&sam.Read1 != 0 {
			fw = fw1
		} else {
			fw = fw2
		}

		if err := fw.WriteRecord(rec.Name, rec.AuxFields, rec.Seq.Expand(), rec.Qual); err != nil {
			return fmt.Errorf("Error writing record %v:\n %v\n", rec.Name, err)
		}
	}
	return nil
}
