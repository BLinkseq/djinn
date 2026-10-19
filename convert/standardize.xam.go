package convert

import "github.com/blinkseq/djinn/xam"

func Standardize(infile string, threads int, asSam bool) error {
	infile, err := xam.FileOrStdin(infile)
	if err != nil {
		return err
	}

	readThread, writeThread := xam.IoThreads(threads)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, r := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := r.Header()
	pg := xam.NewPG(hdr, "djinn standardize "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	// ── open writer ───────────────────────────────────────────────────────────
	outChan, doneChan := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop record channel ──────────────────────────────────────────────

	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal != "" {
			xam.SetBxString(rec, &bxVal)
			xam.SetVX(rec, vxVal)
		}
		// push updated record into writer channel
		outChan <- rec
	}

	close(outChan)

	<-doneChan
	return nil
}
