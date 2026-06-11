package main

import (
	"flag"
	"fmt"
	"os"

	"djinn/xam"
)

func main() {
	threads := flag.Int("threads", 1, "Number of threads for reading BAM. Diminishing returns beyond 2-4.")
	invalid := flag.Bool("invalid", false, "Include invalid barcodes")

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: djinn sam extract [options] input.bam > output.bc\n")
	}
	flag.Parse()
	args := flag.Args()
	infile := xam.FileOrStdin(args, flag.Usage)
	set := make(map[string]struct{}, 5_000_000)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, *threads)

	// ── loop record channel ──────────────────────────────────────────────
	for rec := range recChan {
		bxVal, hasBX, vxVal := xam.FindBarcode(rec)
		if !hasBX {
			continue
		}
		if !vxVal && !*invalid {
			continue
		}
		set[bxVal] = struct{}{}
	}
	for key := range set {
		fmt.Printf("%s\n", key)
	}
}
