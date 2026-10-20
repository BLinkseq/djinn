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
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		invalid, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}
		return count.CountXam(args[0], invalid, threads)
	},
}

// TODO: Docstring is all wrong
var countFqCmd = &cobra.Command{
	Use:     "count [options] file.R1.fq.gz [file.R2.fq.gz]",
	Short:   "Count barcode occurance",
	Example: "count -t 4 curratum.R1.fq.gz curratum.R2.fq.gz > curratum.bc",
	Long: "Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped). " +
		"Inputs expect barcodes to follow the standard (BX tag), stlfr (@seq_id#barcode), or tellseq " +
		"(@seq_id:barcode) formats. Writes to stdout.",
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
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		invalid, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}
		return count.CountFQ(args, invalid)
	},
}

func init() {
	samCmd.AddCommand(countXamCmd)
	fqCmd.AddCommand(countFqCmd)

	//---Command line arguments-------------
	countXamCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
	countXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")

	countFqCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes")
}
