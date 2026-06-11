package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"

	"djinn/xam"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func main() {
	nThreads := flag.Int("threads", 1, "Number of threads for BAM io. Diminishing returns beyond 4.")
	asSam := flag.Bool("sam", false, "Output as uncompressed SAM")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "usage: djinn sam standardize [options] [bam|sam]\n\n")
		fmt.Fprintf(os.Stderr, "options:\n")
		flag.PrintDefaults()
	}
	flag.Parse()

	args := flag.Args()
	infile := xam.FileOrStdin(args, flag.Usage)

	threads := min(max(*nThreads, 1), runtime.NumCPU())
	readThread := 1
	writeThread := 1

	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}

	r, cleanup, err := xam.OpenXAM(infile, readThread, 2<<20)
	if err != nil {
		log.Fatal(err)
	}
	defer cleanup()

	newHeader, err := xam.AddPG(r.Header(), "djinn sam standardize "+infile)
	if err != nil {
		log.Fatal(err)
	}

	var w xam.AlignmentWriter
	switch *asSam {
	case true:
		sw, err := sam.NewWriter(os.Stdout, newHeader, sam.FlagDecimal)
		if err != nil {
			log.Fatalf("creating SAM writer: %v", err)
		}
		w = &xam.SamWriter{Writer: sw}
	case false:
		w, err = bam.NewWriterLevel(os.Stdout, newHeader, 4, writeThread)
		if err != nil {
			log.Fatalf("creating BAM writer: %v", err)
		}
	}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("reading record: %v", err)
		}
		bxVal, hasBX, VX := xam.FindBarcode(rec)
		if hasBX {
			xam.SetBX(rec, bxVal)
			xam.SetVX(rec, VX)
		}
		if err := w.Write(rec); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing BAM: %v\n", err)
			os.Exit(1)
		}
	}
}
