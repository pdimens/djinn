// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/standardize"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var stdXamCmd = &cobra.Command{
	Use:     "standardize [options] file.bam",
	Short:   "Convert standard linked-read format with BX:Z and VX:i tags",
	Example: "standardize -t 10 curratum.bam > curratum.std.bam",
	Long: "This conversion moves the barcode to the `BX:Z` tag in sam/bam records, maintaining the same barcode type. " +
		"See the documentation for a deeper look into the location and format expectations for different linked-read technologies. " +
		"Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid). Writes to stdout.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide an input file")
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
		if err != nil {
			return err
		}
		return standardize.Standardize(args[0], threads, sam)
	},
}

func init() {
	samCmd.AddCommand(stdXamCmd)

	//---Command line arguments-------------
	stdXamCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM")
	stdXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
}
