// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/deconvolve"

	"github.com/spf13/cobra"
)

var deconvolveXamCmd = &cobra.Command{
	Use:     "deconvolve <-d> DISTANCE [options] file.bam",
	Short:   "Tell apart molecules that share a barcode",
	Example: "deconvolve -d 50000 -@ 4 curratum.bam > curratum.deconv.bam",
	Long: "Inputs must be one coordinate-sorted SAM/BAM file, or - to read it from stdin. Writes every record to stdout. " +
		"Where the reads of one barcode are more than DISTANCE bases apart on a chromosome, they are taken to come from " +
		"different molecules, and the barcode (BX tag) of the later ones gets a hyphenated number: BARCODE, then " +
		"BARCODE-1, BARCODE-2, and so on, starting over on each chromosome. Nothing else about the records changes. " +
		"A summary table (molecules before and after, and how many barcodes were split how many times) is written to stderr.",
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
		if args[0] != "-" {
			if err := filecheck(args[0]); err != nil {
				return err
			}
			if err := isSamSorted(args[0]); err != nil {
				return err
			}
		}
		distance, err := cmd.Flags().GetInt("distance-threshold")
		if err != nil {
			return err
		}
		if distance <= 0 {
			return fmt.Errorf("--distance-threshold is required and must be greater than 0, got %d", distance)
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		distance, err := cmd.Flags().GetInt("distance-threshold")
		if err != nil {
			return err
		}
		asSam, err := cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}
		return deconvolve.DeconvolveXam(args[0], distance, asSam, threads)
	},
}

func init() {
	samCmd.AddCommand(deconvolveXamCmd)

	//---Command line arguments-------------
	deconvolveXamCmd.Flags().IntP("distance-threshold", "d", 50000, "Distance (bp) between a barcode's reads that makes them separate molecules")
	deconvolveXamCmd.Flags().BoolP("sam", "S", false, "Output sam instead of BAM")
	deconvolveXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")
}
