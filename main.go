package djinn

import (
	fastqreader "djinn/fastqreader"
	"io"
	"iter"
	"log"
	"os"
	"os/exec"
)

// first, a sentinel function that reads the first 200 records of R1/R2 and looks to assess:
// 1. if it's already standardized
// 2. if not, what the fastq format type is (haplotagging, stlfr, tellseq)

func main(r1, r2 string) {
	var r1_out = "PREFIX.R1.fq.gz"
	var r2_out = "PREFIX.R2.fq.gz"
	var stdinR1 io.WriteCloser
	var stdinR2 io.WriteCloser
	var format string
	var err error
	var barcodes map[string]string

	format, err = fastqreader.FindFastqFormat(r1, r2)
	if err != nil {
		log.Printf("Error opening %s and/or %s to identify file formatting. If the files are unable to be opened at this preprocessing stage, they will be unable to be read for alignment.", r1, r2)
		os.Exit(1)
	}

	// Create the gzip command that will write to output.gz
	cmd_R1 := exec.Command("gzip", "-c")
	cmd_R2 := exec.Command("gzip", "-c")

	// Get the stdin pipe to write data to gzip
	stdinR1, err = cmd_R1.StdinPipe()
	if err != nil {
		log.Fatal("Error creating stdin pipe:", err)
	}
	stdinR2, err = cmd_R2.StdinPipe()
	if err != nil {
		log.Fatal("Error creating stdin pipe:", err)
	}

	// Create the output files
	outFileR1, err := os.Create(r1_out)
	if err != nil {
		log.Fatal("Error creating output file:", err)
	}
	outFileR2, err := os.Create(r2_out)
	if err != nil {
		log.Fatal("Error creating output file:", err)
	}

	defer outFileR1.Close()
	defer outFileR2.Close()

	// Set gzip's stdout to write to the files
	cmd_R1.Stdout = outFileR1
	cmd_R2.Stdout = outFileR2

	// Start the gzip process
	if err := cmd_R1.Start(); err != nil {
		log.Fatal("Error starting gzip: ", err)
	}
	if err := cmd_R2.Start(); err != nil {
		log.Fatal("Error starting gzip: ", err)
	}

	// Important: close stdin when done
	defer stdinR1.Close()
	defer stdinR2.Close()

	var record fastqreader.FastQRecord
	var recordNew fastqreader.FastQRecord

	// no need to check for error b/c file was already opened by sentinel
	fqr, _ := fastqreader.OpenFastQ(r1, r2)

	// create the barcode lookup table
	barcodes = make(map[string]string)

	// create the 20nt barcode generator
	generator := BarcodeGenerator()
	next, stop := iter.Pull(generator)
	defer stop()

	// READ TILL THE END
	for {
		err := fqr.ReadOneLine(&record)
		if err != nil {
			break // I THINK THIS MIGHT BE THE END OF THE FILE?
		}

		bc, ok := barcodes[record.Barcode]
		// if barcode wasn't previously seen...
		if !ok {
			combo, ok := next()
			if !ok {
				log.Fatal("No more unique barcodes available for conversion. Your input data has more than 4^20 unique barcodes.")
			}
			barcodes[record.Barcode] = combo
		}
		// replace the barcode in the FASTQ record to the 20nt version
		record.Barcode = barcodes[record.Barcode]

		recordNew = convertFunc(record)
		// WRITE R1

		// WRITE R2

	}
}
