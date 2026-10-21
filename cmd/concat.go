// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/concat"

	"github.com/spf13/cobra"
)

// preCmd represents the preprocess command
var concatXamCmd = &cobra.Command{
	Use:                   "concat [options] file1.bam file2.bam...",
	Short:                 "Barcode-aware file concatenation",
	Example:               "concat -@ 4 curratum_gulf.bam curratum_gulf2.bam > curratum.concat.bam",
	Long:                  "Concatenate records from linked-read SAM/BAM files while making sure barcodes tags (BX:Z) remain unique for every sample by overwriting all barcodes. Input must be one or more SAM/BAM files. Writes to stdout.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
		}
		if err := cobra.MinimumNArgs(1)(cmd, args); err != nil {
			return fmt.Errorf("Concatenation requires at least 1 input SAM-format file (SAM/BAM)")
		}
		if err := filecheck(args[0]); err != nil {
			return err
		}
		if len(args) == 2 {
			if err := filecheck(args[1]); err != nil {
				return err
			}
		}
		if _, err := cmd.Flags().GetBool("sam"); err != nil {
			return err
		}

		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		asSam, _ := cmd.Flags().GetBool("sam")
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		if len(args) == 1 {
			return concat.ConcatOneXam(args[0], threads, asSam)
		}
		return concat.ConcatXam(args, threads, asSam)
	},
}

func init() {
	samCmd.AddCommand(concatXamCmd)
	//fqCmd.AddCommand(concatFqCmd)

	//---Command line arguments-------------
	concatXamCmd.Flags().BoolP("sam", "S", false, "Output SAM instead of BAM")
	concatXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")

	//sortFqCmd.Flags().IntP("threads", "@", 2, "Threads to use")
	//sortFqCmd.Flags().StringP("tmp-prefix", "t", "", "Folder for temporary files")
	//sortFqCmd.Flags().Bool("standardize", false, "Output using standardized BX:Z:/VX:i: tags instead of the detected input format")
	//sortFqCmd.Flags().StringP("bc-list", "b", "", "Valid barcode list (10X Genomics input format only)")
}
