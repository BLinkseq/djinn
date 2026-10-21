// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/count"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var countXamCmd = &cobra.Command{
	Use:                   "count [options] file.bam",
	Short:                 "Count barcode occurance",
	Example:               "count -@ 4 curratum.bam > curratum.bc",
	Long:                  "Inputs must be one SAM/BAM file. Writes to stdout.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
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
		invalid, _ := cmd.Flags().GetBool("invalid")
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		return count.CountXam(args[0], invalid, threads)
	},
}

var countFqCmd = &cobra.Command{
	Use:                   "count [options] file.R1.fq.gz [file.R2.fq.gz]",
	Short:                 "Count unique barcode occurance",
	Example:               "count curratum.R1.fq.gz curratum.R2.fq.gz > curratum.bc",
	Long:                  "Accepts one FASTQ file or file pair (R1 and R2, can be gzipped). Writes to stdout.",
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
		if err := cobra.MaximumNArgs(2)(cmd, args); err != nil {
			return err
		}
		for i := range args {
			if err := filecheck(args[i]); err != nil {
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
			if err := filecheck(mapfile); err != nil {
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

		return count.CountFQ(args, invalid, mapfile)
	},
}

func init() {
	samCmd.AddCommand(countXamCmd)
	fqCmd.AddCommand(countFqCmd)

	//---Command line arguments-------------
	countXamCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
	countXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")

	countFqCmd.Flags().StringP("bc-list", "b", "", "Valid barcode list [if 10X input]")
	countFqCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
	countFqCmd.Flags().IntP("threads", "@", 2, "Threads to use")
}
