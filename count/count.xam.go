package count

import (
	"bufio"
	"os"
	"strconv"

	"github.com/blinkseq/djinn/xam"
)

func CountXam(infile string, invalid bool, threads int) error {
	infile, err := xam.FileOrStdin(infile)
	if err != nil {
		return err
	}
	set := make(map[string]int, 7_000_000)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, threads)

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
		set[bxVal]++
	}

	for key, val := range set {
		writer.WriteString(key)
		writer.WriteByte('\t')
		writer.WriteString(strconv.Itoa(int(val)))
		writer.WriteByte('\n')
	}
	return nil
}
