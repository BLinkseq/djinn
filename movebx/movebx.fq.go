package movebx

import (
	"bufio"
	"fmt"
	"io"
	"strconv"

	"github.com/blinkseq/djinn/fastq"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func MoveBxFq(fqs []string, prefix string) error {
	var converter func(*fastq.CoreFq, *bufio.Writer)

	coreParser, fmtType, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fqs[0])
	}

	switch fmtType {
	case "standard":
		converter, err = fastq.FormatConverter(fmtType, true)
		if err != nil {
			return err
		}
	case "haplotagging":
		converter, err = fastq.FormatConverter(fmtType, false)
		if err != nil {
			return err
		}
	default:
		return fmt.Errorf("%s linked-read format does not have BX tags to move to the end", fmtType)
	}

	// ── open reader ───────────────────────────────────────────────────────────
	for idx, i := range fqs { // iterate over files
		// ---- FQ reader -------------------------
		fqReader, err := fastx.NewReader(seq.DNA, i, "")
		if err != nil {
			return fmt.Errorf("opening %s: %w", i, err)
		}

		// ---- FQ writer -------------------------
		outfq, err := xopen.Wopen(prefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz")
		if err != nil {
			return err
		}

		// iterate through records
		for {
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			core, _ := coreParser(rec)
			converter(&core, outfq.Writer)
		}
		fqReader.Close()
		outfq.Close()
	}
	return nil
}
