// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/invalid"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

var invalidXamCmd = &cobra.Command{
	Use:     "filter-invalid [options] file.bam",
	Short:   "Remove reads with invalid barcodes",
	Example: "filter-invalid -t 10 curratum.bam > curratum.valid.bam",
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
		inv, err := cmd.Flags().GetString("invalid")
		if err != nil {
			return err
		}
		err = ensureWritableDir(inv)
		if err != nil {
			return err
		}
		return invalid.FilterInvalidXam(args[0], inv, sam, threads)
	},
}

var invalidFqCmd = &cobra.Command{
	Use:     "filter-invalid <-i> PREFIX file.R1.fq file.R2.fq",
	Short:   "Remove reads with invalid barcodes",
	Example: "filter-invalid -i invalid/omykiss omykiss.R1.fq omykiss.R2.fq",
	Long: "Inputs must be one or two FASTQ files (R1 and R2, can be gzipped). " +
		"Inputs expect barcodes to follow the standard (BX tag), stlfr (@seq_id#barcode), or tellseq " +
		"(@seq_id:barcode) formats.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
		}
		if err := cobra.MinimumNArgs(2)(cmd, args); err != nil {
			return err
		}
		if err := cobra.MaximumNArgs(3)(cmd, args); err != nil {
			return err
		}
		inv, err := cmd.Flags().GetString("invalid")
		if err != nil {
			return err
		}
		if inv == args[0] {
			return fmt.Errorf("File prefix for invalid and valid output cannot be the same.")
		}
		for _, i := range args[1:] {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		inv, err := cmd.Flags().GetString("invalid")
		if err != nil {
			return err
		}
		err = ensureWritableDir(args[0])
		if err != nil {
			return err
		}
		err = ensureWritableDir(inv)
		if err != nil {
			return err
		}
		return invalid.FilterInvalidFQ(args[1:], args[0], inv)
	},
}

func init() {
	samCmd.AddCommand(invalidXamCmd)
	fqCmd.AddCommand(invalidFqCmd)

	//---Command line arguments-------------
	invalidXamCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM")
	invalidXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	invalidXamCmd.Flags().StringP("invalid", "i", "", "Write records with invalid barcodes to files with this prefix")

	invalidFqCmd.Flags().StringP("invalid", "i", "", "Write records with invalid barcodes to files with this prefix")
}
