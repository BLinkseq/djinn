package invalid

import (
	"fmt"
	"github.com/blinkseq/djinn/fastq"
	"io"
	"strconv"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func FilterInvalidFQ(fqs []string, prefix, invalidprefix string) error {
	keepInvalid := invalidprefix != ""
	seq.ValidateSeq = false

	// determine what kind of linked-read tech it is, once, from the first file
	coreParser, _, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fqs[0])
	}

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

		var outfqInv *xopen.Writer
		if keepInvalid {
			outfqInv, err = xopen.Wopen(invalidprefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz")
			if err != nil {
				return err
			}
		}

		for { // iterate through records
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

			core, ok := coreParser(rec)
			if !ok || !core.Valid {
				if keepInvalid {
					core.ToStandard(outfqInv.Writer)
				}
				continue
			}

			core.ToStandard(outfq.Writer)
		}
		fqReader.Close()
		outfq.Close()
		if keepInvalid {
			outfqInv.Close()
		}
	}

	return nil
}
