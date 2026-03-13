package extract

import (
	"bufio"
	"flag"
	"fmt"
	"os"

	"djinn/bam/xam"

	"github.com/biogo/hts/bam"
)

func Extract() error {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: djinn bam extract input.bam > output.bc\n")
	}
	args := flag.Args()
	if len(args) != 1 {
		flag.Usage()
		os.Exit(1)
	}
	infile := args[0]

	set := make(map[string]struct{}, 1_000_000)
	// ── open BAM ──────────────────────────────────────────────────────────────
	bf, err := os.Open(infile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening BAM: %v\n", err)
		os.Exit(1)
	}
	defer bf.Close()

	br, err := bam.NewReader(bufio.NewReaderSize(bf, 2<<20), 1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating BAM reader: %v\n", err)
		os.Exit(1)
	}
	defer br.Close()

	// ── loop through BAM records ───────────────────────────────────────────────
	for {
		rec, err := br.Read()
		if err != nil {
			break // EOF
		}

		bxVal, hasBX := xam.GetStringTag(rec, "BX")
		if !hasBX {
			continue
		}
		set[bxVal] = struct{}{}

	}
	for key := range set {
		fmt.Printf("%s\n", key)
	}
	return nil
}
