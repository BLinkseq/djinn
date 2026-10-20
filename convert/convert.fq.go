package convert

import (
	"bufio"
	"fmt"
	"github.com/blinkseq/djinn/barcodes"
	"github.com/blinkseq/djinn/fastq"
	"io"
	"os"
	"slices"
	"strconv"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func ConvertFq(fqs []string, convTo, prefix, bcmap string) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var err error
	var bcs barcodes.Generator
	var converter func(*fastq.CoreFq, *bufio.Writer)
	var fmtType string

	// --- assess the data type and get the parser -----------------------------
	coreParser, fmtType, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		if bcmap == "" {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s. If this is 10X data, a barcode file must be provided to identify inline barcodes.", fqs[0])
		} else {
			// SWAP TO THE 10X CONVERTER and return
			return ConvertFrom10X(fqs, convTo, prefix, bcmap)
		}
	}

	if fmtType == convTo {
		return fmt.Errorf("trying to convert to and from identical format (%s)", convTo)
	}

	// --- create the barcode conversion generator -------------------------------------

	switch convTo {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
		converter = (*fastq.CoreFq).ToHaplotagging
	case "stlfr":
		bcs = barcodes.NewStlfr()
		converter = (*fastq.CoreFq).ToStlfr
	case "tellseq":
		bcs = barcodes.NewTellseq()
		converter = (*fastq.CoreFq).ToTellseq
	case "10x":
		bcs, err = barcodes.NewTenXList()
		if err != nil {
			return err
		}
		converter = (*fastq.CoreFq).ToTenX
	default:
		return fmt.Errorf("unknown barcode type %q", convTo)
	}
	defer bcs.Close()

	seen := make(map[string][]byte, 4_000_000)
	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer

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
			_rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			rec, ok := coreParser(_rec)
			// if no barcode was found, add invalid BX, put to write and continue
			// VX is already invalid when no barcode found, doesn't need to be explicitly declared
			if !ok {
				rec.Barcode = bcs.GetInvalid()
				converter(&rec, outfq.Writer)
				continue
			}

			if convertedBC, isPresent := seen[string(rec.Barcode)]; isPresent {
				// barcode previously seen, pull the converted barcode
				rec.Barcode = convertedBC
			} else {
				// barcode not yet observed, generate new barcode, add it to map
				n, ok := bcs.NextInto(bcBuf)
				if !ok {
					return fmt.Errorf("too many unique barcodes for the conversion technology requested — unable to generate more barcodes.")
				}
				// independent copy; bcBuf gets reused/overwritten next iteration
				newBC := append([]byte(nil), bcBuf[:n]...)
				seen[string(rec.Barcode)] = newBC
				rec.Barcode = newBC
			}
			// write to output buffer
			converter(&rec, outfq.Writer)
		}
		fqReader.Close()
		outfq.Close()
	}

	// ------ write barcodes to map ------------------------
	f, err := os.Create(prefix + ".bc.map")
	if err != nil {
		return err
	}
	defer f.Close()

	writer := bufio.NewWriter(f)
	defer writer.Flush()

	for key, val := range seen {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.Write(val)
		writer.WriteByte('\n')
	}

	return nil
}

// ConvertFrom10X converts 10X-format FASTQ (barcode inline in R1's
// sequence prefix) to another linked-read technology. fqs is either a
// single R1-only file or an [R1, R2] pair -- 10X's barcode lives solely
// in R1, so R2 (when present) only ever needs that barcode propagated
// onto it by matching read ID, never re-derived itself.
func ConvertFrom10X(fqs []string, convTo, prefix, bcmap string) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var bcs barcodes.Generator
	var converter func(*fastq.CoreFq, *bufio.Writer)

	// --- create the barcode conversion generator -------------------------------------
	switch convTo {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
		converter = (*fastq.CoreFq).ToHaplotagging
	case "stlfr":
		bcs = barcodes.NewStlfr()
		converter = (*fastq.CoreFq).ToStlfr
	case "tellseq":
		bcs = barcodes.NewTellseq()
		converter = (*fastq.CoreFq).ToTellseq
	case "10x":
		return fmt.Errorf("trying to convert to and from identical format (10X)")
	}
	defer bcs.Close()

	seen := make(map[string][]byte, 4_000_000)
	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer
	bclist, err := barcodes.LoadTenXFromFile(bcmap)
	if err != nil {
		return err
	}

	// resolveBarcode maps core's inline-extracted (or absent) 10X barcode
	// onto its converted-technology equivalent, generating a fresh one
	// the first time a barcode is seen.
	resolveBarcode := func(core *fastq.CoreFq, ok bool) error {
		if !ok {
			core.Barcode = bcs.GetInvalid()
			return nil
		}
		if convertedBC, isPresent := seen[string(core.Barcode)]; isPresent {
			core.Barcode = convertedBC
			return nil
		}
		n, ok2 := bcs.NextInto(bcBuf)
		if !ok2 {
			return fmt.Errorf("too many unique barcodes for the conversion technology requested — unable to generate more barcodes.")
		}
		// independent copy; bcBuf gets reused/overwritten next iteration
		newBC := append([]byte(nil), bcBuf[:n]...)
		seen[string(core.Barcode)] = newBC
		core.Barcode = newBC
		return nil
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
			if err := resolveBarcode(&r1core, ok); err != nil {
				return err
			}
			converter(&r1core, R1Writer.Writer)
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
				if err := resolveBarcode(&r1core, ok); err != nil {
					return err
				}
				converter(&r1core, R1Writer.Writer)
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
					r2core.Valid = r1core.Valid
				}
				converter(&r2core, R2Writer.Writer)
			default:
				return r2err
			}
		}
	}

	// ------ write barcodes to map ------------------------
	f, err := os.Create(prefix + ".bc.map")
	if err != nil {
		return err
	}
	defer f.Close()

	writer := bufio.NewWriter(f)
	defer writer.Flush()

	for key, val := range seen {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.Write(val)
		writer.WriteByte('\n')
	}

	return nil
}
