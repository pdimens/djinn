// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/count"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var countXamCmd = &cobra.Command{
	Use:                   "count [options] file.bam",
	Short:                 "Count barcode occurance",
	Example:               "count -t 4 curratum.bam > curratum.bc",
	Long:                  "Inputs must be one or two FASTQ files (R1 and R2, can be gzipped). Writes to stdout.",
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

		invalid, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}
		return count.CountXam(args[0], invalid, threads)
	},
}

var countFqCmd = &cobra.Command{
	Use:     "count [options] file.bam",
	Short:   "Count barcode occurance",
	Example: "count -t 4 curratum.bam > curratum.bc",
	Long: "Inputs must one or two FASTQ files (R1 and R2, can be gzipped).Inputs expect barcodes to follow the standard (BX tag), stlfr (@seq_id#barcode), or tellseq " +
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
		for i := range args {
			if err := filecheck(args[i]); err != nil {
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
		return count.CountFQ(args, invalid)
	},
}

func init() {
	samCmd.AddCommand(countXamCmd)
	fqCmd.AddCommand(countFqCmd)

	//---Command line arguments-------------
	countXamCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
	countFqCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
}
