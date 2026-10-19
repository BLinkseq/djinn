// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/stats"
	"fmt"
	"runtime"

	"github.com/spf13/cobra"
)

var statsXamCmd = &cobra.Command{
	Use:     "stats [options] file.bam",
	Short:   "Linked-read metrics per molecule",
	Example: "stats -d 50000 -@ 4 curratum.bam > curratum.molecules.tsv",
	Long: "Inputs must be one coordinate-sorted SAM/BAM file. Writes a table to stdout with one line per molecule: " +
		"number of reads, start, end, inferred length, aligned bp, inferred insert length, and molecule coverage " +
		"(based on aligned bases and on inferred inserts). Barcodes come from the BX tag (or a TELLseq/stLFR-style " +
		"read name), and a barcode is valid if the VX tag says so, or, when there is no VX tag, if the barcode itself " +
		"is valid. Reads without a valid barcode are pooled into one \"invalid\" line per contig. " +
		"Use --distance-threshold to split reads sharing a barcode into separate molecules (BARCODE-1, BARCODE-2, ...) " +
		"when the gap between consecutive reads exceeds it.",
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
		distance, err := cmd.Flags().GetInt("distance-threshold")
		if err != nil {
			return err
		}
		if distance < 0 {
			return fmt.Errorf("--distance-threshold must be 0 or greater, got %d", distance)
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
		return stats.StatsXam(args[0], distance, threads)
	},
}

func init() {
	samCmd.AddCommand(statsXamCmd)

	//---Command line arguments-------------
	statsXamCmd.Flags().IntP("distance-threshold", "d", 0, "Distance (bp) that splits a barcode's reads into separate molecules, 0 = never split")
	statsXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
}
