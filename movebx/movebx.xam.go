package movebx

import (
	"os"
	"strings"

	"github.com/biogo/hts/sam"
	"github.com/blinkseq/djinn/xam"
)

var BXAUX [2]byte = [2]byte{'B', 'X'}

// MoveTagLast shifts the first occurrence of tag (such as "BX") to the end of
// rec.AuxFields in place without allocating a new slice.
func moveTagLast(rec *sam.Record) {
	aux := rec.AuxFields
	for i, entry := range aux {
		if entry.Tag() == BXAUX {
			target := aux[i]
			copy(aux[i:], aux[i+1:])
			aux[len(aux)-1] = target
		}
	}
}

func MoveBxXam(infile string, threads int, asSam bool) error {
	readThread, writeThread := xam.IoThreads(threads)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)
	defer br.Close()

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn "+strings.Join(os.Args[1:], " "))
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop record channel ──────────────────────────────────────────────
	for rec := range recChan {
		moveTagLast(rec)
		writeChan <- rec
	}
	close(writeChan)
	<-writeDone

	return nil
}
