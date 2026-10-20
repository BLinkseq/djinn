package convert

import (
	"fmt"
	"github.com/blinkseq/djinn/barcodes"
	"github.com/blinkseq/djinn/fastq"
	"io"
	"slices"
	"strconv"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func StandardizeFq(fqs []string, prefix, bcmap string) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var err error
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
			return StandardizeFrom10X(fqs, prefix, bcmap)
		}
	}

	// --- pull the invalid BC from the format detected -------------------------------------
	var invalid []byte

	switch fmtType {
	case "haplotagging":
		invalid = barcodes.NewHaplotagging().GetInvalid()
	case "stlfr":
		invalid = barcodes.NewStlfr().GetInvalid()
	case "tellseq":
		invalid = barcodes.NewTellseq().GetInvalid()
	default:
		return fmt.Errorf("unknown format type %q", fmtType)
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
			_rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			rec, ok := coreParser(_rec)
			// if no barcode was found, add invalid BX
			// VX is already invalid when no barcode found, doesn't need to be explicitly declared
			if !ok {
				rec.Barcode = invalid
			}
			rec.ToStandard(outfq.Writer)
		}
		fqReader.Close()
		outfq.Close()
	}
	return nil
}

// StandardizeFrom10X standardizes 10X-format FASTQ (barcode inline in
// R1's sequence prefix). fqs is either a single R1-only file or an
// [R1, R2] pair -- 10X's barcode lives solely in R1, so R2 (when
// present) only ever needs that barcode propagated onto it by matching
// read ID, never re-derived itself.
func StandardizeFrom10X(fqs []string, prefix, bcmap string) error {
	invalid := barcodes.NewTenX().GetInvalid()

	bclist, err := barcodes.LoadTenXFromFile(bcmap)
	if err != nil {
		return err
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
			if !ok {
				r1core.Barcode = invalid
			}
			r1core.ToStandard(R1Writer.Writer)
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

	var r1core fastq.CoreFq
	var r2core fastq.CoreFq

	// iterate through records
Loop:
	for {
		// process R1, which doesn't need an R2
		_r1, r1err := R1Reader.Read()
		switch r1err {
		case io.EOF:

		case nil:
			var ok bool
			r1core, ok = fastq.TenXR12Corefq(_r1, bclist)
			// if no barcode was found, add invalid BX
			// VX is already invalid when no barcode found, doesn't need to be explicitly declared
			if !ok {
				r1core.Barcode = invalid
			}
			r1core.ToStandard(R1Writer.Writer)
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
			r2core.ToStandard(R2Writer.Writer)
		default:
			return r2err
		}
	}
	return nil
}
