package main

import (
	"djinn/bam/xam"
	"flag"
	"fmt"
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
	fmt.Println(xam.DetermineBarcode("TAGAGAN"))
	//args := flag.Args()
	//extract.Extract()
}
