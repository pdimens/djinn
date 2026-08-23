// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/convert"
	"fmt"
	"runtime"
	"strings"

	"github.com/spf13/cobra"
)

var convertXamCmd = &cobra.Command{
	Use:     "convert [options] toFormat file.bam",
	Short:   "Convert barcodes between linked-read formats",
	Example: "convert -t 3 10x gcurratum.bam > gcurratum.10x.bam",
	Long: `Inputs must be one SAM/BAM file, with barcodes specified in BX:Z SAM tag. Writes to stdout. The first positional argument toFormat specifies the target data format (haplotagging, tellseq, stlfr, 10x).

| toFormat     | barcode format | example            |
------------------------------------------------------
| 10x          | 16 nucleotides | GGTTGACATTAAGACA   |
| haplotagging | AxxCxxBxxDxx   | A01C93B56D11       |
| stlfr        | 1_2_3          | 541_9_312          |
| tellseq      | 18 nucleotides | GGCAAATATCGAGAAGTC |`,
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide a target format and input SAM/BAM file")
		}
		if err := cobra.ExactArgs(2)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[1]); err != nil {
			return err
		}
		var err error
		switch strings.ToLower(args[0]) {
		case "haplotagging":
			err = nil
		case "tellseq":
			err = nil
		case "stlfr":
			err = nil
		case "10x":
			err = nil
		default:
			err = fmt.Errorf("toFormat must be one of: haplotagging, stlfr, tellseq, 10x")
		}
		if err != nil {
			return err
		}
		mapfile, err := cmd.Flags().GetString("map")
		if err != nil {
			return err
		}
		if mapfile == "" {
			return fmt.Errorf("conversion map file must be provided via --map/-m")
		}
		if err := filecheck(mapfile); err != nil {
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

		sam, err := cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}
		return convert.ConvertXam(args[1], strings.ToLower(args[0]), threads, sam)
	},
}

func init() {
	samCmd.AddCommand(convertXamCmd)
	//fqCmd.AddCommand(convertFqCmd)

	//---Command line arguments-------------
	convertXamCmd.Flags().BoolP("sam", "S", false, "Output as SAM instead of BAM")
	convertXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	convertXamCmd.Flags().StringP("map", "m", "", "Write barcode conversion map to this file (required)")

	//convertFqCmd.Flags().StringP("singletons", "s", "", "Write valid singleton records to files with this prefix (optional)")
	//convertFqCmd.Flags().StringP("bc-count", "b", "", "Write valid barcode counts to this file (optional)")

}
