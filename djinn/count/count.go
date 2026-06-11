package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"

	"djinn/xam"
)

func main() {
	threads := flag.Int("threads", 1, "Number of threads for BAM io. Diminishing returns beyond 2-4.")
	invalid := flag.Bool("invalid", false, "Include invalid barcodes")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: djinn sam count [options] input.bam > output.bc\n")
	}
	flag.Parse()
	args := flag.Args()
	infile := xam.FileOrStdin(args, flag.Usage)

	set := make(map[string]int16, 5_000_000)
	// ── open BAM ──────────────────────────────────────────────────────────────
	r, cleanup, err := xam.OpenXAM(infile, *threads, 2<<20)
	if err != nil {
		log.Fatal(err)
	}
	defer cleanup()

	// ── loop through records ───────────────────────────────────────────────
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("reading record: %v", err)
		}

		bxVal, hasBX := xam.FindBarcode(rec)
		if !hasBX {
			continue
		}
		if xam.IsInvalid(bxVal) && !*invalid {
			continue
		}
		set[bxVal]++
	}
	for key, val := range set {
		fmt.Printf("%s\t%d\n", key, val)
	}
}
