package ncbi

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"slices"

	"github.com/blinkseq/djinn/barcodes"
	"github.com/blinkseq/djinn/fastq"

	"github.com/biogo/hts/bgzf"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func NcbiFq(fqs []string, bcmap string, threads int) error {
	var err error

	// --- assess the data type and get the parser -----------------------------
	coreParser, _, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		if bcmap == "" {
			return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s. If this is 10X data, a barcode file must be provided to identify inline barcodes.", fqs[0])
		} else {
			// SWAP TO THE 10X CONVERTER and return
			return NCBIfq10x(fqs, bcmap, threads)
		}
	}

	// ---- BAM writer -------------------------
	// One writer shared across every input file: creating a fresh one
	// per file (as this used to) would discard whatever the previous
	// file's writer had buffered but not yet flushed.
	bgzfW := bgzf.NewWriter(os.Stdout, threads) // concurrent blocks
	bw := bufio.NewWriter(bgzfW)
	bw.WriteString("@HD	VN:1.6	SO:unsorted	GO:query\n")

	// ── open reader ───────────────────────────────────────────────────────────
	// iterate over files
	for _, i := range fqs {
		// ---- FQ reader -------------------------
		fqReader, err := fastx.NewReader(seq.DNA, i, "")
		if err != nil {
			return fmt.Errorf("opening %s: %w", i, err)
		}

		// ---- iterate through records -------------------------
		for {
			_rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			rec, _ := coreParser(_rec)
			rec.WriteBam(bw)
		}
		fqReader.Close()
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	if err := bgzfW.Close(); err != nil {
		return err
	}
	if err := os.Stdout.Close(); err != nil {
		return err
	}
	return nil
}

// NCBIfq10x converts 10X-format FASTQ (barcode inline in R1's sequence
// prefix) to unmapped BAM. fqs is either a single R1-only file or an
// [R1, R2] pair -- 10X's barcode lives solely in R1, so R2 (when
// present) only ever needs that barcode propagated onto it by matching
// read ID, never re-derived itself.
func NCBIfq10x(fqs []string, bcmap string, threads int) error {
	// ---- read in barcodes -------------------------------
	bclist, err := barcodes.LoadTenXFromFile(bcmap)
	if err != nil {
		return err
	}

	R1Reader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
	}
	defer R1Reader.Close()

	// ---- BAM writer -------------------------
	bgzfW := bgzf.NewWriter(os.Stdout, threads)
	bw := bufio.NewWriter(bgzfW)
	bw.WriteString("@HD	VN:1.6	SO:unsorted	GO:query\n")

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
			r1core, _ := fastq.TenXR12Corefq(_r1, bclist)
			r1core.WriteBam(bw)
		}
		if err := bw.Flush(); err != nil {
			return err
		}
		if err := bgzfW.Close(); err != nil {
			return err
		}
		if err := os.Stdout.Close(); err != nil {
			return err
		}
		return nil
	}

	R2Reader, err := fastx.NewReader(seq.DNA, fqs[1], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[1], err)
	}
	defer R2Reader.Close()

	var r1core fastq.CoreFq
	var r2core fastq.CoreFq

Loop:
	for {
		// process R1, which doesn't need an R2
		_r1, r1err := R1Reader.Read()
		switch r1err {
		case io.EOF:

		case nil:
			r1core, _ = fastq.TenXR12Corefq(_r1, bclist)
			r1core.WriteBam(bw)
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
			r2core.WriteBam(bw)
		default:
			return r2err
		}
	}

	if err := bw.Flush(); err != nil {
		return err
	}
	if err := bgzfW.Close(); err != nil {
		return err
	}
	if err := os.Stdout.Close(); err != nil {
		return err
	}
	return nil
}
