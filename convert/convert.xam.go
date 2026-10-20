package convert

import (
	"bufio"
	"fmt"
	"github.com/blinkseq/djinn/barcodes"
	"github.com/blinkseq/djinn/xam"
	"os"
)

func ConvertXam(infile, convTo, bcMap string, threads int, asSam bool) error {
	// ── init inventory and generator ──────────────────────────────────────────
	var bcs barcodes.Generator
	var err error

	switch convTo {
	case "haplotagging":
		bcs = barcodes.NewHaplotagging()
	case "stlfr":
		bcs = barcodes.NewStlfr()
	case "tellseq":
		bcs = barcodes.NewTellseq()
	case "10x":
		bcs, err = barcodes.NewTenXList()
		if err != nil {
			return err
		}
	default:
		return fmt.Errorf("unknown barcode type %q", convTo)
	}
	defer bcs.Close()

	bcBuf := make([]byte, bcs.MaxLen()) // generated-barcode buffer

	readThread, writeThread := xam.IoThreads(threads)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)
	defer br.Close()

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam convert "+convTo+" "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop record channel ──────────────────────────────────────────────
	var n int
	var ok bool
	seen := make(map[string][]byte, 4_000_000)
	invalidBarcode := bcs.GetInvalid()
	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal == "" || !vxVal {
			xam.SetBxByte(rec, &invalidBarcode)
			writeChan <- rec
			continue
		}

		if converted, ok := seen[bxVal]; ok {
			xam.SetBxByte(rec, &converted)
			writeChan <- rec
			continue
		}

		// not seen — generate the next conversion barcode
		n, ok = bcs.NextInto(bcBuf)
		if !ok {
			return fmt.Errorf("Barcodes exhausted before conversion finished. The data provided has more unique barcodes than %q barcodes can support.", convTo)
		}
		converted := make([]byte, n) // must copy: buf gets reused next iteration
		copy(converted, bcBuf[:n])
		seen[bxVal] = converted

		xam.SetBxByte(rec, &converted)
		writeChan <- rec
	}
	close(writeChan)
	<-writeDone

	// ------ write barcodes to map ------------------------
	f, err := os.Create(bcMap)
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
