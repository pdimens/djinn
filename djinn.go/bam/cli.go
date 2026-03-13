package main

import (
	"flag"
)

func printHelpCmd(_ []string) error {
	flag.Usage()
	return nil
}

//func usage() {
//	intro := `Manipulate linked-read SAM/BAM files.
//
//Usage:
//  djinn sam <command> options... file.bam`
//	fmt.Fprintln(os.Stderr, intro)
//
//	fmt.Fprintln(os.Stderr, "\nCommands:")
//	for _, cmd := range commands {
//		fmt.Fprintf(os.Stderr, "  %-8s %s\n", cmd.Name, cmd.Help)
//	}
//
//	flag.PrintDefaults()
//
//	fmt.Fprintln(os.Stderr)
//	fmt.Fprintf(os.Stderr, "Run `djinn sam <command> -h` to get help for a specific command\n\n")
//}
