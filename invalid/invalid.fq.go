package invalid

import (
	"bufio"
	"fmt"
	"io"
	"slices"
	"strconv"

	"github.com/blinkseq/djinn/fastq"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func FilterInvalidFQ(fqs []string, prefix string, keepInvalid bool, bcmap string) error {
	if bcmap != "" {
		return FilterInvalid10XFQ(fqs, prefix, keepInvalid, bcmap)
	}
	seq.ValidateSeq = false

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

	// ----- iterate over files ----------------
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

		var outfqInv *xopen.Writer
		if keepInvalid {
			outfqInv, err = xopen.Wopen(prefix + ".invalid.R" + strconv.Itoa(idx+1) + ".fq.gz")
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
					writer(&core, outfqInv.Writer)
				}
				continue
			}

			writer(&core, outfq.Writer)
		}
		fqReader.Close()
		outfq.Close()
		if keepInvalid {
			outfqInv.Close()
		}
	}

	return nil
}

func FilterInvalid10XFQ(fqs []string, prefix string, keepInvalid bool, bcmap string) error {
	seq.ValidateSeq = false
	bclist, err := fastq.Load10X(bcmap)
	if err != nil {
		return err
	}

	// ---- create R1 reader and writer
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
				if keepInvalid {
					r1core.ToTenX(R1Writer.Writer)
				}
				continue
			}
			r1core.ToTenX(R1Writer.Writer)
		}
	} else {
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
					if keepInvalid {
						r1core.ToTenX(R1Writer.Writer)
					}
					continue
				}
				r1core.ToTenX(R1Writer.Writer)
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
				// check if R1 exists and match the seq ID
				if r1err == nil {
					if ok := slices.Equal(r1core.ID, r2core.ID); !ok {
						return fmt.Errorf("Misaligned R1 and R2 records. The 10X format only stores the barcode at the beginning of read 1, meaning Read 2 must have the same sequence ID as Read 1 to reliably assign a barcode to Read 2, but they do not match.\nRead 1:\n%s\nRead 2:\n%s", r1core.ID, r2core.ID)
					}
					r2core.Barcode = r1core.Barcode
					if !r1core.Valid {
						if keepInvalid {
							r2core.ToTenX(R2Writer.Writer)
						}
						continue
					}
				}
			default:
				return r2err
			}
		}
	}
	return nil
}
