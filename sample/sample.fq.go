package sample

import (
	"bufio"
	"fmt"
	"github.com/blinkseq/djinn/fastq"
	"io"
	"math"
	"math/rand"
	"os"
	"slices"
	"strconv"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func SampleFq(fqs []string, prefix string, downsample float64, seed int, keepInvalid bool) error {
	// --- Get barcode map --------------------------
	set := make(map[string]struct{}, 4_000_000)
	seq.ValidateSeq = false

	coreParser, _, err := fastq.CoreFqParser(fqs[0])
	if err != nil {
		return fmt.Errorf("%w", err)
	}
	if coreParser == nil {
		return fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fqs[0])
	}

	for _, fq := range fqs {
		fqReader, err := fastx.NewReader(seq.DNA, fq, "")
		if err != nil {
			return fmt.Errorf("opening %s: %w", fq, err)
		}
		for {
			// iterate through records
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

			core, ok := coreParser(rec)
			if !ok {
				continue
			} // no barcode

			if !core.Valid && !keepInvalid {
				continue
			} // is invalid and was not asked to keep invalid

			bc := string(core.Barcode)
			if _, ok := set[bc]; ok {
				continue
			} else {
				set[bc] = struct{}{}
			}
		}
		fqReader.Close()
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
	slices.Sort(keys) // ensures reproducibility later by setting a fixed pre-permutation state

	// downsample by shuffling barcodes -> selecting first n
	var r *rand.Rand
	if seed >= 0 {
		r = rand.New(rand.NewSource(int64(seed)))
	} else {
		r = rand.New(rand.NewSource(time.Now().UTC().UnixNano()))
	}
	r.Shuffle(nBC, func(i, j int) { keys[i], keys[j] = keys[j], keys[i] })

	var nBcToKeep int
	if downsample < 1.0 {
		nBcToKeep = int(math.Round(float64(nBC) * downsample))
	} else {
		nBcToKeep = int(downsample)
	}

	bcFile, err := os.Create(prefix + ".bc")
	if err != nil {
		return err
	}
	bcWriter := bufio.NewWriter(bcFile)
	bcToKeep := make(map[string]struct{}, nBcToKeep)
	for i := range nBcToKeep {
		bcToKeep[keys[i]] = struct{}{}
		bcWriter.WriteString(keys[i])
		bcWriter.WriteByte('\n')
	}
	bcWriter.Flush()
	err = bcFile.Close()
	if err != nil {
		return err
	}

	for idx, fq := range fqs {
		fqReader, err := fastx.NewDefaultReader(fq)
		if err != nil {
			return fmt.Errorf("opening %s: %w", fq, err)
		}

		// ---- FQ writer -------------------------
		outfq, err := xopen.Wopen(prefix + ".R" + strconv.Itoa(idx+1) + ".fq.gz")
		if err != nil {
			return err
		}

		for {
			// iterate through records
			rec, err := fqReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

			core, ok := coreParser(rec)
			if !ok || (!core.Valid && !keepInvalid) {
				continue
			} // no barcode or is invalid and was not asked to keep invalid

			if _, ok := bcToKeep[string(core.Barcode)]; ok {
				core.ToStandard(outfq.Writer)
			}
		}
		fqReader.Close()
		outfq.Close()
	}
	return nil
}
