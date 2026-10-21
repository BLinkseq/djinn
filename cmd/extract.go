// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/extract"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var extractXamCmd = &cobra.Command{
	Use:                   "extract [options] file.bam",
	Short:                 "Extract all unique barcodes",
	Example:               "extract -@ 4 bombus.bam > bombus.bc",
	Long:                  "Inputs must be one SAM/BAM file. Writes to stdout.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide a SAM/BAM input")
		}
		if err := cobra.ExactArgs(1)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[0]); err != nil {
			return err
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		invalid, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}

		return extract.Extract(args[0], invalid, threads)
	},
}

var extractFqCmd = &cobra.Command{
	Use:                   "extract <-i> file.fq",
	Short:                 "Extract all unique barcodes",
	Example:               "extract bombus.R1.fq bombus.R2.fq > bombus.bc",
	Long:                  "Inputs must be any number of FASTQ files (R1 and R2, can be gzipped). Writes to stdout.",
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
		for _, i := range args {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		mapfile, err := cmd.Flags().GetString("bc-list")
		if err != nil {
			return err
		}
		switch mapfile {
		case "", "10x", "10X": // skip check
		default:
			if filecheck(mapfile) != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		invalid, _ := cmd.Flags().GetBool("invalid")
		mapfile, _ := cmd.Flags().GetString("bc-list")
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		return extract.ExtractFQ(args, invalid, mapfile)
	},
}

func init() {
	samCmd.AddCommand(extractXamCmd)
	fqCmd.AddCommand(extractFqCmd)

	//---Command line arguments-------------
	extractXamCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
	extractXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")

	convertFqCmd.Flags().StringP("bc-list", "b", "", "Valid barcode list [if 10X input]")
	extractFqCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
	extractFqCmd.Flags().IntP("threads", "@", 2, "Threads to use")
}
