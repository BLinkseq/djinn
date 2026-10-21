package count

import (
	"bufio"
	"fmt"
	"io"
	"maps"
	"os"
	"strconv"

	"github.com/blinkseq/djinn/fastq"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func CountFQ(fqs []string, invalid bool, bcmap string) error {
	if bcmap != "" {
		return Count10XFQ(fqs, invalid, bcmap)
	}
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
	set := make(map[string]int, 7_000_000) // barcode container

	// R1
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

		core, ok := coreParser(rec)
		if !ok {
			continue
		} // no barcode

		if !core.Valid && !invalid {
			continue
		} // is invalid and was not asked to report invalid

		set[string(core.Barcode)]++
	}
	fqReader.Close()

	// R2
	if len(fqs) > 1 {
		setR2 := make(map[string]int, 500_000) // R2 barcode container

		fqReader, err = fastx.NewDefaultReader(fqs[1])
		if err != nil {
			return fmt.Errorf("opening %s: %w", fqs[1], err)
		}
		defer fqReader.Close()

		for {
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
			// if the barcode already appears in the R1 set, skip
			if _, ok := set[bc]; ok {
				continue
			}
			setR2[bc]++
		}
		maps.Copy(set, setR2)
	}

	for key, val := range set {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.WriteString(strconv.Itoa(int(val)))
		writer.WriteByte('\n')
	}

	return nil
}

func Count10XFQ(fqs []string, invalid bool, bcmap string) error {
	seq.ValidateSeq = false
	bclist, err := fastq.Load10X(bcmap)
	if err != nil {
		return err
	}

	// ── open barcode writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	// ---- main loop -------------------------
	set := make(map[string]int, 7_000_000) // barcode container

	// R1 -- (R2 doesn't have barcode info)
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

		set[string(core.Barcode)]++
	}
	fqReader.Close()

	for key, val := range set {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.WriteString(strconv.Itoa(int(val)))
		writer.WriteByte('\n')
	}

	return nil
}
