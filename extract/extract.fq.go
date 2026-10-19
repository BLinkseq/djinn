package extract

import (
	"bufio"
	"fmt"
	"github.com/blinkseq/djinn/fastq"
	"io"
	"os"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func ExtractFQ(fqs []string, invalid bool) error {
	// determine what kind of linked-read tech it is
	coreParser, _, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fqs[0])
	}
	seq.ValidateSeq = false
	// ── open barcode writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	// ---- main loop -------------------------
	set := make(map[string]struct{}, 7_000_000) // barcode container

	for _, fq := range fqs {
		if err := extractFromFastq(fq, invalid, coreParser, set, writer); err != nil {
			return err
		}
	}
	return nil
}

// extractFromFastq reads a single FASTQ file and writes newly-seen barcodes
// to writer. It closes its reader before returning, rather than deferring
// the close to the caller's lifetime (which, across many input files, would
// keep every prior file handle open until ExtractFQ itself returned).
func extractFromFastq(fq string, invalid bool, coreParser func(*fastx.Record) (fastq.CoreFq, bool), set map[string]struct{}, writer *bufio.Writer) error {
	fqReader, err := fastx.NewReader(seq.DNA, fq, "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fq, err)
	}
	defer fqReader.Close()

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
	return nil
}
