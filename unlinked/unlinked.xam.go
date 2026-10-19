package unlinked

import (
	"bufio"
	"github.com/blinkseq/djinn/xam"
	"io"
	"os"
	"strconv"

	"github.com/biogo/hts/sam"
)

func getCount(infile string, threads int) map[string]int {
	set := make(map[string]int, 7_000_000)
	// ── open reader ───────────────────────────────────────────────────────────
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, threads)

	// ── loop record channel ──────────────────────────────────────────────
	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal != "" && vxVal {
			set[bxVal]++
		}
	}
	return set
}

func RemoveUnlinkedXam(infile, singletons, barcodecount string, asSam bool, threads int) error {
	// guard against draining stdin when getting barcode counts
	var f *os.File
	if infile == "-" {
		var err error
		f, err = os.CreateTemp("", "xam-spool-*.bam")
		if err != nil {
			return err
		}
		defer os.Remove(f.Name())
		if _, err := io.Copy(f, os.Stdin); err != nil {
			return err
		}
		if _, err := f.Seek(0, io.SeekStart); err != nil {
			return err
		}
		infile = f.Name()
	}

	bcCounts := getCount(infile, threads)
	// write barcode counts if requested
	if barcodecount != "" {
		f, err := os.Create(barcodecount)
		if err != nil {
			return err
		}
		writer := bufio.NewWriter(f)
		for key, val := range bcCounts {
			writer.WriteString(key)
			writer.WriteByte('\t')
			writer.WriteString(strconv.Itoa(int(val)))
			writer.WriteByte('\n')
		}
		writer.Flush()
		f.Close()
	}

	readThread, writeThread := xam.IoThreads(threads)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam filter-singletons "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)
	var singletonChan chan *sam.Record
	var singleDone chan bool

	if singletons != "" {
		singletonChan, singleDone = xam.NewXamWriterChan(singletons, hdr, xam.ChanCap, xam.IoBuf, 1, asSam)
	} else {
		singleDone = make(chan bool)
		close(singleDone)
	}

	// ── loop record channel ──────────────────────────────────────────────

	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if (bxVal == "") || !vxVal {
			continue
		}
		if bcCounts[bxVal] >= 2 {
			writeChan <- rec
		} else if singletonChan != nil {
			singletonChan <- rec
		}
	}
	close(writeChan)
	<-writeDone

	if singletonChan != nil {
		close(singletonChan)
		<-singleDone
	}
	return nil
}
