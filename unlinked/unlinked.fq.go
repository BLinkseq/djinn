package unlinked

import (
	"bufio"
	"fmt"
	"io"
	"maps"
	"os"
	"slices"
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

// getFqCount10X counts valid 10X barcodes from R1 only -- unlike the other
// linked-read formats, a 10X barcode lives solely in R1's sequence prefix,
// so R2 has nothing to independently count.
func getFqCount10X(fqs []string, bclist map[string]struct{}) (map[string]int, error) {
	seq.ValidateSeq = false
	set := make(map[string]int, 4_000_000)

	fqReader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", fqs[0], err)
	}
	defer fqReader.Close()

	for {
		rec, err := fqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		core, ok := fastq.TenXR12Corefq(rec, bclist)
		if !ok || !core.Valid {
			continue
		}
		set[string(core.Barcode)]++
	}
	return set, nil
}

// RemoveUnlinked10XFQ is the 10X-input counterpart of RemoveUnlinkedFQ: a
// barcode is "linked" if its R1 sequence-prefix matches the whitelist on 2
// or more reads, same threshold as every other format.
func RemoveUnlinked10XFQ(fqs []string, prefix string, keepSingle, barcodecount bool, bcmap string) error {
	bclist, err := fastq.Load10X(bcmap)
	if err != nil {
		return err
	}

	bcCounts, err := getFqCount10X(fqs, bclist)
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

	R1Reader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
	}
	defer R1Reader.Close()

	R1Writer, err := xopen.Wopen(prefix + ".R1.fq.gz")
	if err != nil {
		return err
	}
	defer R1Writer.Close()

	var R1WriterSingle *xopen.Writer
	if keepSingle {
		R1WriterSingle, err = xopen.Wopen(prefix + ".unlinked.R1.fq.gz")
		if err != nil {
			return err
		}
		defer R1WriterSingle.Close()
	}

	route := func(core *fastq.CoreFq, w, wSingle *xopen.Writer) {
		switch {
		case bcCounts[string(core.Barcode)] >= 2:
			core.ToTenX(w.Writer)
		case keepSingle:
			core.ToTenX(wSingle.Writer)
		}
	}

	if len(fqs) == 1 {
		// Single-end 10X input: no R2 to propagate the barcode onto.
		for {
			_r1, err := R1Reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			r1core, ok := fastq.TenXR12Corefq(_r1, bclist)
			if !ok || !r1core.Valid {
				continue
			}
			route(&r1core, R1Writer, R1WriterSingle)
		}
		return nil
	}

	R2Reader, err := fastx.NewReader(seq.DNA, fqs[1], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[1], err)
	}
	defer R2Reader.Close()

	R2Writer, err := xopen.Wopen(prefix + ".R2.fq.gz")
	if err != nil {
		return err
	}
	defer R2Writer.Close()

	var R2WriterSingle *xopen.Writer
	if keepSingle {
		R2WriterSingle, err = xopen.Wopen(prefix + ".unlinked.R2.fq.gz")
		if err != nil {
			return err
		}
		defer R2WriterSingle.Close()
	}

	var r1core, r2core fastq.CoreFq

Loop:
	for {
		// process R1, which doesn't need an R2
		_r1, r1err := R1Reader.Read()
		switch r1err {
		case io.EOF:

		case nil:
			var ok bool
			r1core, ok = fastq.TenXR12Corefq(_r1, bclist)
			if !ok || !r1core.Valid {
				r1core.Barcode = nil // marks this pair as having no usable barcode for R2 below
			} else {
				route(&r1core, R1Writer, R1WriterSingle)
			}
		default:
			return r1err
		}

		// process R2
		_r2, r2err := R2Reader.Read()
		switch r2err {
		case io.EOF:
			// break when both are exhausted
			if r1err == io.EOF {
				break Loop
			}
		case nil:
			r2core = fastq.TenXR22Corefq(_r2)
			if r1err == nil && len(r1core.Barcode) > 0 {
				if ok := slices.Equal(r1core.ID, r2core.ID); !ok {
					return fmt.Errorf("Misaligned R1 and R2 records. The 10X format only stores the barcode at the beginning of read 1, meaning Read 2 must have the same sequence ID as Read 1 to reliably assign a barcode to Read 2, but they do not match.\nRead 1:\n%s\nRead 2:\n%s", r1core.ID, r2core.ID)
				}
				r2core.Barcode = r1core.Barcode
				route(&r2core, R2Writer, R2WriterSingle)
			}
		default:
			return r2err
		}
	}
	return nil
}

func RemoveUnlinkedFQ(fqs []string, prefix string, keepSingle, barcodecount bool, bcmap string) error {
	if bcmap != "" {
		return RemoveUnlinked10XFQ(fqs, prefix, keepSingle, barcodecount, bcmap)
	}
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
