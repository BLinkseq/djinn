package unlinked

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
	"github.com/shenwei356/xopen"
)

func getFqCount(infiles []string) (map[string]int, error) {
	set := make(map[string]int, 4_000_000) // barcode container
	var emptymap map[string]int            // for returning empty thing on error
	seq.ValidateSeq = false

	// determine what kind of linked-read tech it is
	coreParser, _, err := fastq.CoreFqParser(infiles[0])
	if err != nil {
		return emptymap, err
	}
	if coreParser == nil {
		return emptymap, fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", infiles[0])
	}

	// R1
	fqReader, err := fastx.NewReader(seq.DNA, infiles[0], "")
	if err != nil {
		return emptymap, fmt.Errorf("opening %s: %w", infiles[0], err)
	}

	for {
		// iterate through records
		rec, err := fqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return emptymap, err
		}
		core, ok := coreParser(rec)
		if !ok || !core.Valid {
			continue
		} // no barcode or is invalid
		set[string(core.Barcode)]++
	}
	fqReader.Close()

	// R2
	if len(infiles) > 1 {
		setR2 := make(map[string]int, 500_000) // R2 barcode container

		fqReader, err = fastx.NewDefaultReader(infiles[1])
		if err != nil {
			return emptymap, fmt.Errorf("opening %s: %w", infiles[1], err)
		}
		defer fqReader.Close()

		for {
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return emptymap, err
			}

			core, ok := coreParser(rec)
			if !ok || !core.Valid {
				continue
			} // no barcode or is invalid
			bc := string(core.Barcode)
			if _, ok := set[bc]; ok {
				continue
			} // if the barcode already appears in the R1 set, skip
			setR2[bc]++
		}
		maps.Copy(set, setR2)
	}
	return set, nil
}

func RemoveUnlinkedFQ(fqs []string, prefix string, keepSingle, barcodecount bool, bcmap string) error {
	// guard against draining stdin when getting barcode counts
	bcCounts, err := getFqCount(fqs)
	if err != nil {
		return err
	}
	if barcodecount {
		f, err := os.Create(prefix + ".count")
		if err != nil {
			return err
		}
		writer := bufio.NewWriter(f)
		for key, val := range bcCounts {
			writer.WriteString(key)
			writer.WriteByte('\t')
			writer.WriteString(strconv.Itoa(val))
			writer.WriteByte('\n')
		}
		writer.Flush()
		f.Close()
	}

	// determine what kind of linked-read tech it is, once, from the first file
	coreParser, bcFmt, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fqs[0])
	}

	// --- determine output format from input ------------
	var writer func(*fastq.CoreFq, *bufio.Writer)

	switch bcFmt {
	case "haplotagging":
		writer = (*fastq.CoreFq).ToHaplotagging
	case "stlfr":
		writer = (*fastq.CoreFq).ToStlfr
	case "tellseq":
		writer = (*fastq.CoreFq).ToTellseq
	}

	// ----- iterate over files --------------------
	for idx, i := range fqs {
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

		var outfqSingle *xopen.Writer
		if keepSingle {
			outfqSingle, err = xopen.Wopen(prefix + ".unlinked.R" + strconv.Itoa(idx+1) + ".fq.gz")
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
				continue
			}
			if bcCounts[string(core.Barcode)] >= 2 {
				writer(&core, outfq.Writer)
			} else if keepSingle {
				writer(&core, outfqSingle.Writer)
			}
		}
		fqReader.Close()
		outfq.Close()
		if keepSingle {
			outfqSingle.Close()
		}
	}

	return nil
}
