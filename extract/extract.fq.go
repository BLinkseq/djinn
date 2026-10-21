package extract

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/blinkseq/djinn/fastq"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func ExtractFQ(fqs []string, invalid bool, bcmap string) error {
	var coreParser func(rec *fastx.Record) (core fastq.CoreFq, ok bool)
	var err error
	if bcmap != "" {
		return Extract10XFQ(fqs, invalid, bcmap)
	} else {
		// determine what kind of linked-read tech it is
		coreParser, _, err = fastq.CoreFqParser(fqs[0])
		if err != nil {
			return fmt.Errorf("%w", err)
		}
		if coreParser == nil {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fqs[0])
		}
	}
	seq.ValidateSeq = false
	// ── open barcode writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	// ---- main loop -------------------------
	set := make(map[string]struct{}, 7_000_000) // barcode container

	for _, fq := range fqs {
		fqReader, err := fastx.NewReader(seq.DNA, fq, "")
		if err != nil {
			return fmt.Errorf("opening %s: %w", fq, err)
		}

		for {
			// iterate through records
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

			core, ok := coreParser(rec)
			if !ok {
				continue
			} // no barcode

			if !core.Valid && !invalid {
				continue
			} // is invalid and was not asked to report invalid

			bc := string(core.Barcode)
			if _, ok := set[bc]; ok {
				continue
			} else {
				set[bc] = struct{}{}
				writer.WriteString(bc)
				writer.WriteByte('\n')
			}
		}
		fqReader.Close()
	}
	return nil
}

func Extract10XFQ(fqs []string, invalid bool, bcmap string) error {
	seq.ValidateSeq = false
	bclist, err := fastq.Load10X(bcmap)
	if err != nil {
		return err
	}

	// ── open barcode writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	set := make(map[string]struct{}, 7_000_000) // barcode container

	fqReader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
	}

	for {
		// iterate through records
		rec, err := fqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		core, ok := fastq.TenXR12Corefq(rec, bclist)
		if !ok {
			continue
		} // no barcode

		if !core.Valid && !invalid {
			continue
		} // is invalid and was not asked to report invalid

		bc := string(core.Barcode)
		if _, ok := set[bc]; ok {
			continue
		} else {
			set[bc] = struct{}{}
			writer.WriteString(bc)
			writer.WriteByte('\n')
		}
	}
	fqReader.Close()
	return nil
}
