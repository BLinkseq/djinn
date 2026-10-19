// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"djinn/stats"
	"fmt"
	"runtime"
	"strings"

	"github.com/spf13/cobra"
)

var statsXamCmd = &cobra.Command{
	Use:     "stats [options] file.bam",
	Short:   "Linked-read metrics per molecule",
	Example: "stats -@ 4 curratum.bam > curratum.molecules.tsv\ndeconvolve -d 50000 curratum.bam | stats - > curratum.molecules.tsv",
	Long: "Inputs must be one coordinate-sorted SAM/BAM file, or - to read it from stdin. Writes a table to stdout with " +
		"one line per molecule (a barcode on a chromosome): number of reads, start, end, inferred length, aligned bp, " +
		"inferred insert length, and molecule coverage (based on aligned bases and on inferred inserts). Barcodes come " +
		"from the BX tag (or a TELLseq/stLFR-style read name), and a barcode is valid if the VX tag says so, or, when " +
		"there is no VX tag, if the barcode itself is valid. Reads without a valid barcode are pooled into one " +
		"\"invalid\" line per contig. To tell apart molecules that share a barcode, run the file through " +
		"`sam deconvolve` first: each numbered barcode (BARCODE-1, BARCODE-2, ...) is then reported as its own molecule.",
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
		if args[0] == "-" {
			return nil // stdin
		}
		return filecheck(args[0])
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))
		return stats.StatsXam(args[0], threads)
	},
}

var statsFqCmd = &cobra.Command{
	Use:     "stats [options] file.fq.gz",
	Short:   "Barcode statistics for a FASTQ file",
	Example: "stats curratum.R1.fq.gz > curratum.R1.stats",
	Long: "Accepts one FASTQ file (can be gzipped). Inputs expect barcodes to follow the standard (BX tag), stlfr " +
		"(@seq_id#barcode), or tellseq (@seq_id:barcode) formats. Writes to stdout the number of reads, the number " +
		"of reads with a barcode, how many of those barcodes are valid and invalid, and then how many barcodes were bad " +
		"at each position: haplotagging segments A, C, B, D (e.g. A00), stlfr parts 1-3 (a 0), or tellseq bases (an N). " +
		"The technology is worked out from the first 100 reads unless --platform is given.",
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
		platform, err := cmd.Flags().GetString("platform")
		if err != nil {
			return err
		}
		switch strings.ToLower(platform) {
		case "", stats.Haplotagging, stats.Stlfr, stats.Tellseq:
			return nil
		}
		return fmt.Errorf("unknown platform %q, expected haplotagging, stlfr or tellseq", platform)
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		platform, err := cmd.Flags().GetString("platform")
		if err != nil {
			return err
		}
		return stats.StatsFQ(args[0], platform)
	},
}

func init() {
	samCmd.AddCommand(statsXamCmd)
	fqCmd.AddCommand(statsFqCmd)

	//---Command line arguments-------------
	statsXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")

	statsFqCmd.Flags().StringP("platform", "p", "", "Linked-read technology: haplotagging, stlfr or tellseq (default: detect from the first 100 reads)")
}
