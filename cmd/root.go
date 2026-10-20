// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"os"

	cc "github.com/ivanpirog/coloredcobra"
	"github.com/spf13/cobra"
)

var VERSION string = "3.0"

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:     "djinn",
	Version: VERSION,
	Short:   "Convert between linked-read formats and barcode styles",
	Long:    "Use the subcommands (e.g. `djinn ...`) to explore the options for FASTQ and SAM files.",
}

var samCmd = &cobra.Command{
	Use:   "sam",
	Short: "SAM/BAM operations",
}

var fqCmd = &cobra.Command{
	Use:   "fq",
	Short: "FASTQ operations",
}

func init() {
	rootCmd.CompletionOptions.DisableDefaultCmd = true
	rootCmd.SetHelpCommand(&cobra.Command{Hidden: true})

	rootCmd.AddCommand(fqCmd)
	rootCmd.AddCommand(samCmd)
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	cc.Init(&cc.Config{
		RootCmd:         rootCmd,
		Headings:        cc.Blue, // + cc.Underline,
		Commands:        cc.HiMagenta + cc.Bold,
		Flags:           cc.HiMagenta, // + cc.Bold,
		FlagsDataType:   cc.Italic,
		NoExtraNewlines: true,
		//Example:         cc.Italic,
		//ExecName:        cc.Bold,
	})
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}
