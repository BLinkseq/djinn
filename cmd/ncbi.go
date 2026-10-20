// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/ncbi"

	"github.com/spf13/cobra"
)

var ncbiXam2FqCmd = &cobra.Command{
	Use:                   "ncbi [options] file.bam",
	Short:                 "BAM → FASTQ conversion from NCBI",
	Example:               "ncbi -@ 10 obesus_stlfr obesus.bam",
	Long:                  "Converts an unmapped SAM/BAM file to FASTQ sequences, losslessly, preserving barcode tags.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide inputs")
		}
		//TODO not exact args, needs min/max
		if err := cobra.ExactArgs(2)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[1]); err != nil {
			return err
		}
		err := ensureWritableDir(args[0])
		if err != nil {
			return err
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))
		return ncbi.NcbiXam(args[1], args[0], threads)
	},
}

var ncbiFq2XamCmd = &cobra.Command{
	Use:                   "ncbi [options] file.R1.fq.gz [file.R2.fq.gz]",
	Short:                 "FASTQ → BAM conversion for NCBI",
	Example:               "ncbi -t 10 obesus.R1.fq.gz obesus.R2.fq.gz > obesus.bam",
	Long:                  "Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped). Converts losslessly to an unmapped SAM/BAM file, preserving barcode tags.",
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
		for _, i := range args {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		mapfile, err := cmd.Flags().GetString("bc-map")
		if err != nil {
			return err
		}
		if mapfile != "" && filecheck(mapfile) != nil {
			return err
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))
		mapfile, _ := cmd.Flags().GetString("bc-map")

		return ncbi.NcbiFq(args, mapfile, threads)
	},
}

// THIS IS PRETTY INCOMPLETE
func init() {
	samCmd.AddCommand(ncbiXam2FqCmd)
	fqCmd.AddCommand(ncbiFq2XamCmd)

	//---Command line arguments-------------
	ncbiXam2FqCmd.Flags().IntP("threads", "@", 2, "Threads to use")

	ncbiFq2XamCmd.Flags().IntP("threads", "@", 2, "Threads to use")
	ncbiFq2XamCmd.Flags().StringP("bc-map", "b", "", "Valid barcode list (10X Genomics input format only)")
}
