// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"os"

	"github.com/spf13/cobra"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:     "djinn",
	Version: "3",
	Short:   "Convert between linked-read formats and barcode styles",
	Long:    "Use the subcommands (e.g. `djinn sam ...`) to explore the options for FASTQ and SAM files.",
}

func init() {
	rootCmd.CompletionOptions.DisableDefaultCmd = true
	rootCmd.SetHelpCommand(&cobra.Command{Hidden: true})
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}
