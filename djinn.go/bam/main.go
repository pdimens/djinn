package main

import (
	"djinn/bam/internal/extract"
	"flag"
	"os"
)

func main() {
	//flag.Usage = usage
	flag.Parse()

	// user needs to provide a subcommand
	if len(flag.Args()) < 1 {
		flag.Usage()
		os.Exit(1)
	}
	//args := flag.Args()
	extract.Extract()
}
