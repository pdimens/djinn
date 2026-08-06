// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/sort"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var sortCmd = &cobra.Command{
	Use:     "sort [options] file.bam",
	Short:   "Sort reads by barcode",
	Example: "sort -t 4 curratum.bam > curratum.sort.bam",
	Long: "Inputs must be one SAM/BAM file or two FASTQ files (R1 and R2, can be gzipped). Both FASTQ and SAM/BAM " +
		"inputs expect barcodes to follow the standard (BX tag), stlfr (@seq_id#barcode), or tellseq " +
		"(@seq_id:barcode) formats. Writes to stdout.",
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
		if len(args) == 2 {
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

		tmpDir, err := cmd.Flags().GetString("tmp-prefix")
		if err != nil {
			return err
		}
		asSam, err := cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}
		return sort.SortByBX(args[0], "-", tmpDir, threads, asSam)
	},
}

func init() {
	rootCmd.AddCommand(sortCmd)

	//---Command line arguments-------------
	sortCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM") //not implemented yet
	sortCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	sortCmd.Flags().StringP("tmp-prefix", "t", "", "Folder for temporary files")
}
