// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/singletons"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var singletonXamCmd = &cobra.Command{
	Use:     "filter-singletons [options] file.bam",
	Short:   "Filter out singleton barcodes",
	Example: "filter-singletons -t 10 curratum.bam > curratum.linked.bam",
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

		sam, err := cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}
		singles, err := cmd.Flags().GetString("singletons")
		if err != nil {
			return err
		}
		err = ensureWritableDir(singles)
		if err != nil {
			return err
		}
		return singletons.FilterSingletonsXam(args[0], singles, sam, threads)
	},
}

func init() {
	samCmd.AddCommand(singletonXamCmd)

	//---Command line arguments-------------
	singletonXamCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM")
	singletonXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	singletonXamCmd.Flags().StringP("singletons", "s", "", "Write valid singleton records to this file")
}
