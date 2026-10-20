package extract

import (
	"bufio"
	"os"

	"github.com/blinkseq/djinn/xam"
)

func Extract(infile string, invalid bool, threads int) error {
	infile, err := xam.FileOrStdin(infile)
	if err != nil {
		return err
	}
	set := make(map[string]struct{}, 4_000_000)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, threads)
	defer br.Close()

	// ── open writer ───────────────────────────────────────────────────────────
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	// ── loop record channel ──────────────────────────────────────────────
	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal == "" {
			continue
		}
		if !vxVal && !invalid {
			continue
		}
		if _, ok := set[bxVal]; ok {
			continue
		} else {
			set[bxVal] = struct{}{}
			writer.WriteString(bxVal)
			writer.WriteByte('\n')
		}
	}
	return nil
}
