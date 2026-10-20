package sample

import (
	"bufio"
	"fmt"
	"github.com/blinkseq/djinn/xam"
	"io"
	"math"
	"math/rand"
	"os"
	"slices"
	"time"
)

func SampleXam(infile string, downsample float64, seed, threads int, invalid, asSam bool) error {
	readThread, writeThread := xam.IoThreads(threads)

	// --- Get barcode map --------------------------
	infile, err := xam.FileOrStdin(infile)
	if err != nil {
		return err
	}

	// guard against draining stdin on the first of two passes over infile
	if infile == "-" {
		f, err := os.CreateTemp("", "xam-spool-*.bam")
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

	set := make(map[string]struct{}, 4_000_000)

	// ── open reader ───────────────────────────────────────────────────────────
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── loop record channel to get barcodes ───────────────────────────────────
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
		}
	}

	// ---- Sort and extract barcodes to keep -----------------------------------
	nBC := len(set)

	if float64(nBC) < downsample {
		return fmt.Errorf("The input has fewer barcodes (%v) than the requested downsampling amount (%v)", nBC, downsample)
	}

	keys := make([]string, 0, nBC)
	for k := range set {
		keys = append(keys, k)
	}
	slices.Sort(keys) // ensures reproducibility later

	// downsample by shuffling barcodes -> selecting first n
	var r *rand.Rand
	if seed >= 0 {
		r = rand.New(rand.NewSource(int64(seed)))
	} else {
		r = rand.New(rand.NewSource(time.Now().UTC().UnixNano()))
	}
	r.Shuffle(nBC, func(i, j int) { keys[i], keys[j] = keys[j], keys[i] })

	var nBcToKeep = int(downsample)
	// write alongside infile, not into the process's cwd: path.Base(infile)
	// dropped infile's directory entirely, writing wherever the command
	// happened to be run from.
	bcFile, err := os.Create(infile + ".bc")
	if err != nil {
		return err
	}

	if downsample < 1.0 {
		nBcToKeep = int(math.Round(float64(nBC) * downsample))
	}
	bcToKeep := make(map[string]struct{}, nBcToKeep)
	bcWriter := bufio.NewWriter(bcFile)
	for i := range nBcToKeep {
		bcToKeep[keys[i]] = struct{}{}
		bcWriter.WriteString(keys[i])
		bcWriter.WriteByte('\n')
	}
	if err := bcWriter.Flush(); err != nil {
		return err
	}

	err = bcFile.Close()
	if err != nil {
		return err
	}

	// ── open reader (again) ───────────────────────────────────────────────────
	recChan, br = xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam sample "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}
	// ── open writer ───────────────────────────────────────────────────────────
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// ── loop record channel to get barcodes ───────────────────────────────────
	for rec := range recChan {
		bxVal, vxVal := xam.FindBarcode(rec)
		if bxVal == "" {
			continue
		}
		if !vxVal && !invalid {
			continue
		}
		if _, ok := bcToKeep[bxVal]; ok {
			writeChan <- rec
		}
	}
	close(writeChan)
	<-writeDone

	return nil
}
