// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"djinn/extract"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var extractCmd = &cobra.Command{
	Use:                   "extract [options] file.bam",
	Short:                 "Extract all unique barcodes",
	Example:               "extract -t 4 bombus.bam > bombus.bc",
	Long:                  "Inputs must be one SAM/BAM file or two FASTQ files (R1 and R2, can be gzipped). Writes to stdout.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
		}
		//TODO not exact args, needs min/max
		if err := cobra.ExactArgs(1)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[0]); err != nil {
			return err
		}
		if len(args) == 3 {
			if err := filecheck(args[1]); err != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, err := cmd.Flags().GetInt("threads")
		if err != nil {
			return err
		}
		maxCores := runtime.NumCPU()
		// clamp between 1 and max system threads
		threads = min(maxCores, max(threads, 1))
		runtime.GOMAXPROCS(threads)

		invalid, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}

		return extract.Extract(args[0], invalid, threads)
	},
}

var extractFqCmd = &cobra.Command{
	Use:     "extract-fq <-i> file.fq",
	Short:   "Extract all unique barcodes",
	Example: "extract bombus.R1.fq bombus.R2.fq > bombus.bc",
	Long: "Inputs must be any number of FASTQ files (R1 and R2, can be gzipped). " +
		"Inputs expect barcodes to follow the standard  (BX tag), stlfr (@seq_id#barcode), or tellseq " +
		"(@seq_id:barcode) formats. Writes to stdout.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
		}

		if err := cobra.MinimumNArgs(1)(cmd, args); err != nil {
			return err
		}
		for _, i := range args {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		invalid, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}
		return extract.ExtractFQ(args, invalid)
	},
}

func init() {
	rootCmd.AddCommand(extractCmd)
	rootCmd.AddCommand(extractFqCmd)

	//---Command line arguments-------------
	extractCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
	extractCmd.Flags().IntP("threads", "@", 2, "Decompression threads to use")

	extractFqCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
}
