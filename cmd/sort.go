// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/sort"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var sortXamCmd = &cobra.Command{
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
		if err := cobra.ExactArgs(1)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[0]); err != nil {
			return err
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

var sortFqCmd = &cobra.Command{
	Use:                   "sort [options] PREFIX file.R1.fq [file.R2.fq]",
	Short:                 "Sort read(s) by barcode",
	Example:               "sort -@ 4 omykiss.sort omykiss.R1.fq.gz omykiss.R2.fq.gz",
	Long:                  "Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped) sharing matching read IDs. Writes PREFIX.R1.fq.gz (and PREFIX.R2.fq.gz if paired).",
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
		if err := ensureWritableDir(args[0]); err != nil {
			return err
		}
		for _, i := range args[1:] {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		bcMap, err := cmd.Flags().GetString("bc-map")
		if err != nil {
			return err
		}
		if bcMap != "" {
			if err := filecheck(bcMap); err != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		tmpDir, _ := cmd.Flags().GetString("tmp-prefix")
		standardize, _ := cmd.Flags().GetBool("standardize")
		bcMap, _ := cmd.Flags().GetString("bc-map")

		return sort.SortFqByBX(args[1:], args[0], tmpDir, threads, standardize, bcMap)
	},
}

func init() {
	samCmd.AddCommand(sortXamCmd)
	fqCmd.AddCommand(sortFqCmd)

	//---Command line arguments-------------
	sortXamCmd.Flags().BoolP("sam", "S", false, "Output SAM instead of BAM") //not implemented yet
	sortXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")
	sortXamCmd.Flags().StringP("tmp-prefix", "t", "", "Folder for temporary files")

	sortFqCmd.Flags().IntP("threads", "@", 2, "Threads to use")
	sortFqCmd.Flags().StringP("tmp-prefix", "t", "", "Folder for temporary files")
	sortFqCmd.Flags().Bool("standardize", false, "Output using standardized BX:Z:/VX:i: tags instead of the detected input format")
	sortFqCmd.Flags().StringP("bc-map", "b", "", "Valid barcode list (10X Genomics input format only)")
}
