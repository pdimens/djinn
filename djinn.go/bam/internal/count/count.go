package count

import (
	"bufio"
	"flag"
	"fmt"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func getStringTag(r *sam.Record, tag string) (string, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			if s, ok := aux.Value().(string); ok {
				return s, true
			}
		}
	}
	return "", false
}

func count(args []string) error {
	flagSet := flag.NewFlagSet("count", flag.ExitOnError)
	flagSet.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: count input.bam > output.bc\n")
		flagSet.PrintDefaults()
	}
	flagSet.Parse()
	args = flagSet.Args()
	if len(args) != 1 {
		flagSet.Usage()
		os.Exit(1)
	}
	infile := args[0]

	set := make(map[string]int16, 1_000_000)
	// ── open BAM ──────────────────────────────────────────────────────────────
	bf, err := os.Open(infile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening BAM: %v\n", err)
		os.Exit(1)
	}
	defer bf.Close()

	br, err := bam.NewReader(bufio.NewReaderSize(bf, 1<<20), 1)
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

		bxVal, hasBX := getStringTag(rec, "BX")
		if !hasBX {
			continue
		}
		set[bxVal]++
	}
	for key, val := range set {
		fmt.Printf("%s\t%d\n", key, val)
	}
	return nil
}
