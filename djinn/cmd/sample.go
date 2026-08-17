// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/sample"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var sampleXamCmd = &cobra.Command{
	Use:     "sample [options] file.bam",
	Short:   "Downsample data by barcode",
	Example: "sample -t 8 -d 50000 curratum.bam > curratum.sample50k.bam",
	Long: "Downsamples a SAM/BAM file by barcode to keep all records containing `-d` randomly sampled barcodes. " +
		"If `d >= 1`, the downsampling is a fixed number of barcodes, whereas `d < 1` would indicate a fraction of the total " +
		"number of barcodes `(e.g. `-d 0.5` retains 50% of all barcodes). Use `--invalid/-i` to include invalid barcodes " +
		"in downsampling. Barcode must be in `BX:Z` SAM tag. Writes to stdout.",
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
		invalid, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}
		downsample, err := cmd.Flags().GetFloat64("downsample")
		if err != nil {
			return err
		}
		if downsample <= 0.0 {
			//cmd.Usage()
			return fmt.Errorf("downsample value must be > 0")
		}
		seed, err := cmd.Flags().GetInt("seed")
		if err != nil {
			return err
		}
		if seed < 0 && seed != -1 {
			fmt.Println("Notice: a negative seed value does not set a random seed, thus is not reproducible")
		}
		if err != nil {
			return err
		}
		return sample.Sample(args[0], downsample, seed, threads, invalid, sam)
	},
}

func init() {
	samCmd.AddCommand(sampleXamCmd)

	//---Command line arguments-------------
	sampleXamCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM")
	sampleXamCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes in subsampling")
	sampleXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	sampleXamCmd.Flags().IntP("seed", "s", -1, "Random seed for sampling, must be >= 0")
	sampleXamCmd.Flags().Float64P("downsample", "d", 0.0, "Number (>=1) or fraction (<1) of barcodes to retain")

}
