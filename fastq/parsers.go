package fastq

import (
	"errors"
	"fmt"
	"io"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// Scan a fastq file and return the type of reading parser it needs, i.e., a haplotagging one,
// tellseq, stlfr, 10x
func CoreFqParser(fq string) (func(rec *fastx.Record) (core CoreFq, ok bool), string, error) {
	var rec *fastx.Record
	var h, t, s int
	var totalReads int
	seq.ValidateSeq = false

	fqReader, err := fastx.NewReader(seq.DNA, fq, "")
	if err != nil {
		return nil, "", fmt.Errorf("opening %s: %w", fq, err)
	}
	defer fqReader.Close()
	for i := range 100 {
		rec, err = fqReader.Read()
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			return nil, "", fmt.Errorf("reading %v, record %v: %w", fq, i, err)
		}
		// Take the paired-end marker (a "/1" or "/2" on the ID, or "1:N:0:ATCG" in the description) out
		// first, as the record parsers do. A barcode ends the ID, so with "/1" still attached
		// (the usual stLFR style: ...#12_34_56/1) it would never be found.
		id, desc, _ := splitCASAVA(rec.ID, rec.Desc)
		// is there a BX tag in the comments/description?
		if _, _, _, _, found := findField(desc, BXTAG); found {
			h += 1
		}
		// if not, look for tellseq
		if _, _, _, found := findTellseqBarcode(id); found {
			t += 1
		}
		// if not, look for stlfr
		if _, _, _, found := findStlfrBarcode(id); found {
			s += 1
		}
		totalReads += 1
	}
	// if more than one style found, return an error
	formatsFound := 0
	if h > 0 {
		formatsFound++
	}
	if s > 0 {
		formatsFound++
	}
	if t > 0 {
		formatsFound++
	}
	if formatsFound > 1 {
		return nil, "", fmt.Errorf("more than one linked-read technology format identified. Input data must use a single format. Reads types identified: Haplotagging - %d | stLFR - %d | TELLseq - %d.", h, s, t)
	}
	switch {
	case h > 0:
		//fmt.Fprintln(os.Stderr, "Format detected: Haplotagging")
		return Haplotag2Corefq, "haplotagging", nil
	case s > 0:
		//fmt.Fprintln(os.Stderr, "Format detected: stLFR")
		return Stlfr2Corefq, "stlfr", nil
	case t > 0:
		//fmt.Fprintln(os.Stderr, "Format detected: TELL-seq")
		return Tellseq2Corefq, "tellseq", nil
	default:
		return nil, "", nil
	}
}
