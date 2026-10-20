// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"fmt"
	"runtime"
	"strings"

	"github.com/blinkseq/djinn/convert"

	"github.com/spf13/cobra"
)

var convertXamCmd = &cobra.Command{
	Use:     "convert [options] toFormat file.bam",
	Short:   "Convert barcodes between linked-read formats",
	Example: "convert -@ 3 -m gcur.map 10x gcurratum.bam > gcurratum.10x.bam",
	Long: `Input must be one SAM/BAM file, with barcodes specified in BX:Z SAM tag. Writes to stdout. The first positional argument toFormat specifies the target data format (standard, haplotagging, tellseq, stlfr, 10x).

| toFormat     | barcode format | example            |
------------------------------------------------------
| 10x          | 16 nucleotides | GGTTGACATTAAGACA   |
| haplotagging | AxxCxxBxxDxx   | A01C93B56D11       |
| stlfr        | 1_2_3          | 541_9_312          |
| tellseq      | 18 nucleotides | GGCAAATATCGAGAAGTC |
| standard     | same as input  | GGCAAATATCGAGAAGTC |
`,
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide a target format and input SAM/BAM file")
		}
		if err := cobra.ExactArgs(2)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[1]); err != nil {
			return err
		}
		var err error
		switch strings.ToLower(args[0]) {
		case "haplotagging", "tellseq", "stlfr", "10x", "standard":
		default:
			return fmt.Errorf("toFormat must be one of: standard, haplotagging, stlfr, tellseq, 10x")
		}
		mapfile, err := cmd.Flags().GetString("map")
		if err != nil {
			return err
		}
		// only raise error for missing map file if conversion not 'standard'
		if strings.ToLower(args[0]) != "standard" {
			if mapfile == "" {
				return fmt.Errorf("conversion map file must be provided via --map/-m")
			}
			if err := ensureWritableDir(mapfile); err != nil {
				return err
			}
		}
		_, err = cmd.Flags().GetInt("threads")
		if err != nil {
			return err
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		runtime.GOMAXPROCS(safethreads(threads))

		mapfile, _ := cmd.Flags().GetString("map")
		sam, err := cmd.Flags().GetBool("sam")
		if err != nil {
			return err
		}

		switch strings.ToLower(args[0]) {
		case "standard":
			return convert.Standardize(args[1], threads, sam)
		default:
			return convert.ConvertXam(args[1], strings.ToLower(args[0]), mapfile, threads, sam)
		}
	},
}

var convertFqCmd = &cobra.Command{
	Use:     "convert [options] outFormat prefix file.R1.fq.gz [file.R2.fq.gz]",
	Short:   "Convert between linked-read FASTQ/barcode formats",
	Example: "convert 10x gcurratum_as10X gcurratum.R1.fq.gz gcurratum.R2.fq.gz",
	Long: `Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped). Writes prefix.R1.fq.gz (and prefix.R2.fq.gz if paired). The first positional argument outFormat specifies the target data format (standard, haplotagging, tellseq, stlfr, 10x). Specifying 10X output will draw barcodes from the official 10X Genomics barcode whitelist.

| outFormat    | barcode format | example            |
------------------------------------------------------
| 10x          | 16 nucleotides | GGTTGACATTAAGACA   |
| haplotagging | AxxCxxBxxDxx   | A01C93B56D11       |
| stlfr        | 1_2_3          | 541_9_312          |
| tellseq      | 18 nucleotides | GGCAAATATCGAGAAGTC |
| standard     | same as input  | GGCAAATATCGAGAAGTC |
`,
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s", cmd.UsageString())
			return fmt.Errorf("please provide a target format and input SAM/BAM file")
		}
		if err := cobra.MinimumNArgs(3)(cmd, args); err != nil {
			return err
		}
		if err := cobra.MaximumNArgs(4)(cmd, args); err != nil {
			return err
		}
		for _, i := range args[2:] {
			if err := filecheck(i); err != nil {
				return err
			}
		}
		var err error
		switch strings.ToLower(args[0]) {
		case "haplotagging", "tellseq", "stlfr", "10x", "standard":
		default:
			return fmt.Errorf("toFormat must be one of: standard, haplotagging, stlfr, tellseq, 10x")
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
		mapfile, _ := cmd.Flags().GetString("bc-map")

		switch strings.ToLower(args[0]) {
		case "standard":
			return convert.StandardizeFq(args[2:], args[1], mapfile)
		default:
			return convert.ConvertFq(
				args[2:],
				strings.ToLower(args[0]),
				args[1],
				mapfile,
			)
		}
	},
}

func init() {
	samCmd.AddCommand(convertXamCmd)
	fqCmd.AddCommand(convertFqCmd)

	//---Command line arguments-------------
	convertXamCmd.Flags().BoolP("sam", "S", false, "Output sam instead of BAM")
	convertXamCmd.Flags().IntP("threads", "@", 2, "Threads to use")
	convertXamCmd.Flags().StringP("map", "m", "", "Write barcode conversion map to this file (required if not 'standard')")

	convertFqCmd.Flags().StringP("bc-map", "b", "", "Valid barcode list (10X Genomics input format only)")
}

//TODO SET THREADS FOR FASTQ
