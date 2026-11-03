package sample

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"slices"
	"strconv"
	"time"

	"github.com/blinkseq/djinn/fastq"

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

func SampleFq10x(fqs []string, prefix string, bcMap string, downsample float64, seed int) error {
	// --- Get barcode map --------------------------
	set := make(map[string]struct{}, 4_000_000)
	seq.ValidateSeq = false
	bclist, err := fastq.Load10X(bcMap)
	if err != nil {
		return err
	}

	// ----- Get barcodes from the R1 records ----------------
	fqReader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
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

		core, ok := fastq.TenXR12Corefq(rec, bclist)
		if !ok {
			continue
		} // no barcode

		if !core.Valid {
			continue
		}
		bc := string(core.Barcode)
		if _, ok := set[bc]; ok {
			continue
		} else {
			set[bc] = struct{}{}
		}
	}
	fqReader.Close()

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
	switch {
	case downsample < 1.0:
		nBcToKeep = int(math.Round(float64(nBC) * downsample))
	default:
		nBcToKeep = int(downsample)
	}
	bcToKeep := make(map[string]struct{}, nBcToKeep)

	defer func() error {
		bcFile, err := os.Create(prefix + ".bc")
		if err != nil {
			return err
		}
		bcWriter := bufio.NewWriter(bcFile)
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
		return nil
	}()

	// ------ Process the R1 and (possibly) R2 together
	R1Reader, err := fastx.NewReader(seq.DNA, fqs[0], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[0], err)
	}
	defer R1Reader.Close()

	R1Writer, err := xopen.Wopen(prefix + ".R1.fq.gz")
	if err != nil {
		return err
	}
	defer R1Writer.Close()

	if len(fqs) == 1 {
		// Single-end 10X input: no R2 to propagate the barcode onto.
		for {
			_r1, err := R1Reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			r1core, ok := fastq.TenXR12Corefq(_r1, bclist)
			if !ok {
				continue
			}
			if _, ok := bcToKeep[string(r1core.Barcode)]; ok {
				r1core.ToTenX(R1Writer.Writer)
			}
		}
		return nil
	}

	// --- instantiate R2 reader and writer -----------------
	R2Reader, err := fastx.NewReader(seq.DNA, fqs[1], "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fqs[1], err)
	}
	defer R2Reader.Close()

	R2Writer, err := xopen.Wopen(prefix + ".R2.fq.gz")
	if err != nil {
		return err
	}
	defer R2Writer.Close()

	var r1core, r2core fastq.CoreFq
	var writeRecord bool

Loop:
	for {
		// process R1, which doesn't need an R2
		_r1, r1err := R1Reader.Read()
		switch r1err {
		case io.EOF:
			break Loop // cannot process R2 without R1
		case nil:
			r1core, ok := fastq.TenXR12Corefq(_r1, bclist)
			if !ok {
				continue Loop
			}
			if _, ok := bcToKeep[string(r1core.Barcode)]; ok {
				r1core.ToTenX(R1Writer.Writer)
				writeRecord = true
			} else {
				writeRecord = false
			}
		default:
			return r1err
		}

		// process R2
		_r2, r2err := R2Reader.Read()
		switch r2err {
		case io.EOF:
			break Loop
		case nil:
			r2core = fastq.TenXR22Corefq(_r2)
			if ok := slices.Equal(r1core.ID, r2core.ID); !ok {
				return fmt.Errorf("Misaligned R1 and R2 records. The 10X format only stores the barcode at the beginning of read 1, meaning Read 2 must have the same sequence ID as Read 1 to reliably assign a barcode to Read 2, but they do not match.\nRead 1:\n%s\nRead 2:\n%s", r1core.ID, r2core.ID)
			}
			if writeRecord {
				r2core.ToTenX(R2Writer.Writer)
			}
		default:
			return r2err
		}
	}
	return nil
}
