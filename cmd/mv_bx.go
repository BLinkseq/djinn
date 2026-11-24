// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

var mvBxXamCmd = &cobra.Command{
	Use:                   "mv-bx [options] file.bam",
	Short:                 "Move BX tags to record ends",
	Example:               "mv-bx -@ 5 versicolor.bam > versicolor.bxterm.bam",
	Long:                  "Inputs must be one SAM/BAM file. Writes to stdout.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide an input alignment file")
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

		sam, err := cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}
		return movebx.moveBxXam(args[0], sam, threads)
	},
}

var mvBxFqCmd = &cobra.Command{
	Use:     "rm-invalid <-i> PREFIX file.R1.fq [file.R2.fq]",
	Short:   "Remove reads with invalid barcodes",
	Example: "rm-invalid -i invalid/omykiss omykiss.valid omykiss.R1.fq omykiss.R2.fq",
	Long: "Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped). " +
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
		err := ensureWritableDir(args[0])
		if err != nil {
			return err
		}
		for _, i := range args[1:] {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))
		return nil
	},
}

func init() {
	samCmd.AddCommand(mvBxXamCmd)
	fqCmd.AddCommand(mvBxFqCmd)

	//---Command line arguments-------------
	invalidXamCmd.Flags().BoolP("sam", "S", false, "Output sam instead of BAM")
	invalidXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")

	//invalidFqCmd.Flags().String("10x", "", "Input is 10X format")
	//invalidFqCmd.Flags().Lookup("10x").NoOptDefVal = "10x" // value used when bare --10x
	//invalidFqCmd.Flags().BoolP("invalid", "i", false, "Write records with invalid barcodes to prefix.invalid")
	//invalidFqCmd.Flags().IntP("threads", "@", 2, "Threads to use")
}

//moveBxXam(infile string, threads int, asSam bool)
