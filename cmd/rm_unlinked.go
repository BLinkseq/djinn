// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"

	"github.com/blinkseq/djinn/unlinked"

	"github.com/spf13/cobra"
)

var unlinkedXamCmd = &cobra.Command{
	Use:     "rm-unlinked [options] file.bam",
	Short:   "Remove records with unlinked barcodes",
	Example: "rm-unlinked -t 10 curratum.bam > curratum.linked.bam",
	Long: "Input must be one SAM/BAM file " +
		"Input expect barcodes to follow the standard (BX tag), stlfr (seq_id#barcode), or tellseq " +
		"(seq_id:barcode) formats. By default, only writes filtered output. Writes to stdout.",
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
		singles, err := cmd.Flags().GetString("singletons")
		if err != nil {
			return err
		}
		if singles != "" {
			if err := ensureWritableDir(singles); err != nil {
				return err
			}
		}
		_, err = cmd.Flags().GetInt("threads")
		if err != nil {
			return err
		}
		_, err = cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}

		bcCount, err := cmd.Flags().GetString("bc-count")
		if err != nil {
			return err
		}
		if bcCount != "" {
			err = ensureWritableDir(bcCount)
			if err != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		sam, _ := cmd.Flags().GetBool("sam")
		singles, _ := cmd.Flags().GetString("singletons")
		bcCount, _ := cmd.Flags().GetString("bc-count")

		return unlinked.RemoveUnlinkedXam(args[0], singles, bcCount, sam, threads)
	},
}

var unlinkedFqCmd = &cobra.Command{
	Use:     "rm-unlinked [options] PREFIX file.R1.fq [file.R2.fq]",
	Short:   "Remove reads with unlinked barcodes",
	Example: "rm-unlinked rclamitans.linked rclamitans.R1.fq rclamitans.R2.fq",
	Long: "Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped). " +
		"Inputs expect barcodes to follow the standard (BX tag), stlfr (@seq_id#barcode), or tellseq " +
		"(@seq_id:barcode) formats. By default, only writes filtered output.",
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
		singles, err := cmd.Flags().GetString("singletons")
		if err != nil {
			return err
		}
		if singles == args[0] {
			return fmt.Errorf("File prefix for singletons and valid output cannot be the same.")
		}
		if singles != "" {
			err = ensureWritableDir(singles)
			if err != nil {
				return err
			}
		}
		err = ensureWritableDir(args[0])
		if err != nil {
			return err
		}
		for _, i := range args[1:] {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		bcCount, err := cmd.Flags().GetString("bc-count")
		if err != nil {
			return err
		}
		if bcCount != "" {
			err = ensureWritableDir(bcCount)
			if err != nil {
				return err
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		singles, _ := cmd.Flags().GetString("singletons")
		bcCount, _ := cmd.Flags().GetString("bc-count")

		return unlinked.RemoveUnlinkedFQ(args[1:], args[0], singles, bcCount)
	},
}

func init() {
	samCmd.AddCommand(unlinkedXamCmd)
	fqCmd.AddCommand(unlinkedFqCmd)

	//---Command line arguments-------------
	unlinkedXamCmd.Flags().BoolP("sam", "S", false, "Output sam instead of BAM")
	unlinkedXamCmd.Flags().IntP("threads", "@", 2, "Worker threads to use")
	unlinkedXamCmd.Flags().StringP("singletons", "s", "", "Write valid singleton records to this file (optional)")
	unlinkedXamCmd.Flags().StringP("bc-count", "b", "", "Write valid barcode counts to this file (optional)")

	unlinkedFqCmd.Flags().StringP("singletons", "s", "", "Write valid singleton records to files with this prefix (optional)")
	unlinkedFqCmd.Flags().StringP("bc-count", "b", "", "Write valid barcode counts to this file (optional)")

}
