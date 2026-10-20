package invalid

import (
	"github.com/blinkseq/djinn/xam"

	"github.com/biogo/hts/sam"
)

func FilterInvalidXam(infile, invalid string, asSam bool, threads int) error {
	readThread, writeThread := xam.IoThreads(threads)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam filter-invalid "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)
	var invalidChan chan *sam.Record
	var invalidDone chan bool

	if invalid != "" {
		invalidChan, invalidDone = xam.NewXamWriterChan(invalid, hdr, xam.ChanCap, xam.IoBuf, 1, asSam)
	} else {
		invalidDone = make(chan bool)
		close(invalidDone)
	}

	// ── loop record channel ──────────────────────────────────────────────
	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if (bxVal != "") && vxVal {
			writeChan <- rec
		} else if invalidChan != nil {
			invalidChan <- rec
		}
	}
	close(writeChan)
	<-writeDone

	if invalidChan != nil {
		close(invalidChan)
		<-invalidDone
	}
	return nil
}
