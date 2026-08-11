// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/ncbi"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var ncbiCmd = &cobra.Command{
	Use:                   "ncbi [options] file.bam",
	Short:                 "BAM → FASTQ conversion from NCBI",
	Example:               "ncbi -t 10 obesus_stlfr obesus.bam",
	Long:                  "Converts an unmapped SAM/BAM file to FASTQ sequences, losslessly, preserving barcode tags.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
		}
		//TODO not exact args, needs min/max
		if err := cobra.ExactArgs(2)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[1]); err != nil {
			return err
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
		//TODO ensure prefix is writable
		//err = ensureWritableDir(cmd.ar)
		if err != nil {
			return err
		}
		return ncbi.NCBI(args[1], args[0], threads)
	},
}

func init() {
	rootCmd.AddCommand(ncbiCmd)

	//---Command line arguments-------------
	ncbiCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM")
	ncbiCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	ncbiCmd.Flags().StringP("singletons", "s", "", "Write valid singleton records to this file")
}
