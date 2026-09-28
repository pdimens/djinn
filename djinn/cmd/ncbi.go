// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/ncbi"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var ncbiXam2FqCmd = &cobra.Command{
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
		err := ensureWritableDir(args[0])
		if err != nil {
			return err
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))
		return ncbi.NcbiXam(args[1], args[0], threads)
	},
}

// THIS IS PRETTY INCOMPLETE
func init() {
	samCmd.AddCommand(ncbiXam2FqCmd)

	//---Command line arguments-------------
	ncbiXam2FqCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM")
	ncbiXam2FqCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	ncbiXam2FqCmd.Flags().StringP("singletons", "s", "", "Write valid singleton records to this file")
}
