package ncbi

import (
	"fmt"
	"github.com/blinkseq/djinn/fastq"
	"github.com/blinkseq/djinn/xam"

	"github.com/biogo/hts/sam"
)

func NcbiXam(infile, prefix string, threads int) error {
	r1Path := prefix + ".R1.fq.gz"
	r2Path := prefix + ".R2.fq.gz"

	readThread, writeThread := xam.IoThreads(threads)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)
	defer br.Close()

	// ── FASTQ writers ──────────────────────────────────────────────
	threadsPerWriter := max(writeThread/2, 1)
	fw1, err := fastq.NewFastqWriter(r1Path, "/1", 4, threadsPerWriter)
	if err != nil {
		return fmt.Errorf("Error opening R1 output: %v\n", err)
	}

	fw2, err := fastq.NewFastqWriter(r2Path, "/2", 4, threadsPerWriter)
	if err != nil {
		fw1.Close()
		return fmt.Errorf("Error opening R2 output: %v\n", err)
	}

	// ── loop record channel ──────────────────────────────────────────────
	var fw *fastq.FastqWriter

	for rec := range recChan {
		if rec.Flags&sam.Read1 != 0 {
			fw = fw1
		} else {
			fw = fw2
		}

		if err := fw.WriteRecord(rec.Name, rec.AuxFields, rec.Seq.Expand(), rec.Qual); err != nil {
			fw1.Close()
			fw2.Close()
			return fmt.Errorf("Error writing record %v:\n %v\n", rec.Name, err)
		}
	}
	if err := fw1.Close(); err != nil {
		return fmt.Errorf("Error closing R1 output: %v\n", err)
	}
	if err := fw2.Close(); err != nil {
		return fmt.Errorf("Error closing R2 output: %v\n", err)
	}
	return nil
}
