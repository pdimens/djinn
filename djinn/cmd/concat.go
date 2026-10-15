// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/sort"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var concatXamCmd = &cobra.Command{
	Use:                   "sort [options] file.bam",
	Short:                 "Sort reads by barcode",
	Example:               "sort -t 4 curratum.bam > curratum.sort.bam",
	Long:                  "Input must be one SAM/BAM file. Writes to stdout.",
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
		_, err := cmd.Flags().GetString("tmp-prefix")
		if err != nil {
			return err
		}
		_, err = cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}

		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		tmpDir, _ := cmd.Flags().GetString("tmp-prefix")
		asSam, _ := cmd.Flags().GetBool("sam")

		return sort.SortByBX(args[0], "-", tmpDir, threads, asSam)
	},
}

func init() {
	samCmd.AddCommand(concatXamCmd)
	//fqCmd.AddCommand(concatFqCmd)

	//---Command line arguments-------------
	concatXamCmd.Flags().BoolP("sam", "S", false, "Output SAM instead of BAM") //not implemented yet
	concatXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")

	//sortFqCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	//sortFqCmd.Flags().StringP("tmp-prefix", "t", "", "Folder for temporary files")
	//sortFqCmd.Flags().Bool("standardize", false, "Output using standardized BX:Z:/VX:i: tags instead of the detected input format")
	//sortFqCmd.Flags().StringP("bc-map", "b", "", "Valid barcode list (10X Genomics input format only)")
}
