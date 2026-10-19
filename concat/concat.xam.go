package concat

import (
	"fmt"
	"github.com/blinkseq/djinn/barcodes"
	"github.com/blinkseq/djinn/xam"
	"os"
	"strconv"
	"strings"
)

func ConcatXam(infiles []string, threads int, asSam bool) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var bcs barcodes.Generator
	var err error

	bcFmt, err := xam.DetectBarcodeFormat(infiles[0])
	if err != nil {
		return err
	}

	switch bcFmt {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
	case "stlfr":
		bcs = barcodes.NewStlfr()
	case "nucleotide-18":
		bcs = barcodes.NewTellseq()
	case "nucleotide-16":
		bcs = barcodes.NewTenX()
	default:
		// "nucleotide-<length>", length is not always 2 digits
		_, lenStr, ok := strings.Cut(bcFmt, "-")
		if !ok {
			return fmt.Errorf("unrecognized barcode format: %s", bcFmt)
		}
		n, err := strconv.Atoi(lenStr)
		if err != nil {
			return fmt.Errorf("unrecognized barcode format %q: %w", bcFmt, err)
		}
		bcs, err = barcodes.NewGeneric(n)
		if err != nil {
			return err
		}
	}
	invalidBarcode := bcs.GetInvalid()
	defer bcs.Close()

	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer

	readThread, writeThread := xam.IoThreads(threads)

	// ── open first file to get header ────────────────────────────────────────────
	_, br := xam.NewXamReaderChan(infiles[0], xam.ChanCap, xam.IoBuf, readThread)

	// update PG line in header
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam concat "+strings.Join(os.Args[1:], " "))
	if err := hdr.AddProgram(pg); err != nil {
		br.Close()
		return err
	}
	br.Close()

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop through remaining files ──────────────────────────────────────────────
	var n int
	var _ok bool

	for _, samfile := range infiles {
		recChan, br := xam.NewXamReaderChan(samfile, xam.ChanCap, xam.IoBuf, readThread)
		seenLocal := make(map[string][]byte, 2_000_000)

		for rec := range recChan {
			bxVal, vxVal := xam.FindBarcode(rec)
			if bxVal == "" {
				// no barcode, set invalid and write
				xam.SetBxByte(rec, &invalidBarcode)
				writeChan <- rec
				continue
			}
			if !vxVal {
				// invalid, just write
				writeChan <- rec
				continue
			}

			// is it seen in the local inventory?
			if _, ok := seenLocal[bxVal]; ok {
				// previously seen, set to the new barcode it was previously assigned
				_bc := seenLocal[bxVal]
				xam.SetBxByte(rec, &_bc)
			} else {
				// doesn't yet exist in local inventory, add with new unique barcode
				n, _ok = bcs.NextInto(bcBuf)
				if !_ok {
					return fmt.Errorf("The data provided has more unique barcodes than the input technology (%q) can support. Stopped at file %q.", bcFmt, samfile)
				}
				converted := make([]byte, n) // must copy: buf gets reused next iteration
				copy(converted, bcBuf[:n])
				seenLocal[bxVal] = converted
				xam.SetBxByte(rec, &converted)
			}
			writeChan <- rec
		}
		br.Close()
	}
	close(writeChan)
	<-writeDone
	return nil
}

// Fake concatenation on a single file (mostly a no-op) by updating the header and writing unmodified records to stdout.
func ConcatOneXam(infile string, threads int, asSam bool) error {
	readThread, writeThread := xam.IoThreads(threads)

	// ── open reader --------------────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)
	defer br.Close()

	// update PG line in header
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam concat "+strings.Join(os.Args[1:], " "))
	if err := hdr.AddProgram(pg); err != nil {
		br.Close()
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	for rec := range recChan {
		writeChan <- rec
	}
	close(writeChan)
	<-writeDone
	return nil
}
