// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/sample"

	"github.com/spf13/cobra"
)

var sampleXamCmd = &cobra.Command{
	Use:     "sample [options] file.bam",
	Short:   "Downsample data by barcode",
	Example: "sample -t 8 -d 50000 curratum.bam > curratum.sample50k.bam",
	Long: "Downsamples a SAM/BAM file by barcode to keep all records containing `-d` randomly sampled barcodes. " +
		"If `d >= 1`, the downsampling is a fixed number of barcodes, whereas `d < 1` would indicate a fraction of the total " +
		"number of barcodes (e.g. `-d 0.5` retains 50% of all barcodes). Use `--invalid/-i` to include invalid barcodes " +
		"in downsampling. Barcode must be in `BX:Z` SAM tag. Writes to stdout.",
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
		_, err := cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}
		_, err = cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}
		downsample, err := cmd.Flags().GetFloat64("downsample")
		if err != nil {
			return err
		}
		if downsample <= 0.0 {
			return fmt.Errorf("downsample value must be > 0")
		}
		_, err = cmd.Flags().GetInt("seed")
		if err != nil {
			return err
		}

		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		sam, _ := cmd.Flags().GetBool("sam")
		invalid, _ := cmd.Flags().GetBool("invalid")
		downsample, _ := cmd.Flags().GetFloat64("downsample")
		seed, _ := cmd.Flags().GetInt("seed")
		if seed < 0 && seed != -1 {
			fmt.Println("Notice: a negative seed value does not set a random seed, thus is not reproducible")
		}
		return sample.SampleXam(args[0], downsample, seed, threads, invalid, sam)
	},
}

var sampleFqCmd = &cobra.Command{
	Use:     "sample [options] prefix file.R1.fq [file.R2.fq]",
	Short:   "Downsample data by barcode",
	Example: "sample -d 50000 gbiloba.50k gbiloba.R1.fq gbiloba.R2.fq",
	Long: "Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped). " +
		"Downsamples by barcode to keep all records containing `-d` randomly sampled barcodes. " +
		"If `d >= 1`, the downsampling is a fixed number of barcodes, whereas `d < 1` would indicate a fraction of the total " +
		"number of barcodes (e.g. `-d 0.5` retains 50% of all barcodes). Use `--invalid/-i` to include invalid barcodes " +
		"in downsampling.",
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
		for _, i := range args[1:] {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		_, err := cmd.Flags().GetBool("invalid")
		if err != nil {
			return err
		}
		downsample, err := cmd.Flags().GetFloat64("downsample")
		if err != nil {
			return err
		}
		if downsample <= 0.0 {
			//cmd.Usage()
			return fmt.Errorf("downsample value must be > 0")
		}
		_, err = cmd.Flags().GetInt("seed")
		if err != nil {
			return err
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		invalid, _ := cmd.Flags().GetBool("invalid")
		downsample, _ := cmd.Flags().GetFloat64("downsample")
		seed, _ := cmd.Flags().GetInt("seed")
		if seed < 0 && seed != -1 {
			fmt.Println("Notice: a negative seed value does not set a random seed, thus is not reproducible")
		}
		return sample.SampleFq(args[1:], args[0], downsample, seed, invalid)
	},
}

func init() {
	samCmd.AddCommand(sampleXamCmd)
	fqCmd.AddCommand(sampleFqCmd)

	//---Command line arguments-------------
	sampleXamCmd.Flags().BoolP("sam", "S", false, "Output sam instead of BAM")
	sampleXamCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes in subsampling")
	sampleXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")
	sampleXamCmd.Flags().IntP("seed", "s", -1, "Random seed for sampling, must be >= 0")
	sampleXamCmd.Flags().Float64P("downsample", "d", 0.0, "Number (>=1) or fraction (<1) of barcodes to keep")

	sampleFqCmd.Flags().BoolP("invalid", "i", false, "Include invalid barcodes in subsampling")
	sampleFqCmd.Flags().IntP("seed", "s", -1, "Random seed for sampling, must be >= 0")
	sampleFqCmd.Flags().Float64P("downsample", "d", 0.0, "Number (>=1) or fraction (<1) of barcodes to keep")

}
