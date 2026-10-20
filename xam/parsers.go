package xam

import (
	"bufio"
	"fmt"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// hasHyphenSuffix reports whether rest is empty (no suffix) or begins
// with '-' followed by at least one more character. Every barcode shape
// matcher in this file (and isStlfrBody/findInlineBarcode in sam.go)
// uses this to permit an optional trailing tag after a technology's own
// barcode shape -- e.g. Cell Ranger's "-1" GEM well suffix convention,
// generalized here to every linked-read technology: a hyphen, then any
// non-whitespace content, is always allowed after the recognized core.
func hasHyphenSuffix(rest string) bool {
	return len(rest) == 0 || (rest[0] == '-' && len(rest) > 1)
}

// isHaplotaggingBarcode reports whether bc has the fixed haplotagging
// shape "A##C##B##D##" -- 12 bytes: a letter followed by two digits,
// repeated for A, C, B, D in that order -- optionally followed by "-"
// and a suffix (see hasHyphenSuffix).
func isHaplotaggingBarcode(bc string) bool {
	if len(bc) < 12 {
		return false
	}
	letters := [4]byte{'A', 'C', 'B', 'D'}
	for i, want := range letters {
		pos := i * 3
		if bc[pos] != want || !isDigit(bc[pos+1]) || !isDigit(bc[pos+2]) {
			return false
		}
	}
	return hasHyphenSuffix(bc[12:])
}

// nucleotideCoreLen reports the length of bc's leading run of DNA bases
// (including the ambiguity code N) -- the shape shared by both TELLseq
// and 10X barcodes once standardized into a BX:Z tag -- when that run
// covers the whole of bc or is followed by a hyphen suffix (see
// hasHyphenSuffix); 0 otherwise. The two technologies are told apart
// purely by this core length (10X: 16, TELLseq: 18), not bc's raw
// length, so a "-1"-style suffix never skews the count; see
// DetectBarcodeFormat.
func nucleotideCoreLen(bc string) int {
	i := 0
	for i < len(bc) && isBase(bc[i]) {
		i++
	}
	if i == 0 || !hasHyphenSuffix(bc[i:]) {
		return 0
	}
	return i
}

// isNucleotideBarcode reports whether bc is a nucleotide barcode, same
// shape nucleotideCoreLen matches.
func isNucleotideBarcode(bc string) bool {
	return nucleotideCoreLen(bc) > 0
}

// openAlignmentReader opens file (or stdin for "-") as a SAM or BAM
// stream, auto-detecting which from the first two bytes. It's meant for
// a short, throwaway read (see DetectBarcodeFormat): SAM/BAM streams
// can't be rewound, so a caller that goes on to actually process the
// file must open a fresh reader of its own afterward -- this one should
// just be closed when done sniffing.
//
// Closing isBAM's fh is the caller's responsibility when isBAM is true:
// *bam.Reader's own Close() only releases its internal decompression
// resources, not the underlying file handle (matching the asymmetry
// NewXamReaderChan's goroutine already works around the same way; SAM's
// wrapper Close() does close its file handle, via SamReader.closer).
func openAlignmentReader(file string) (r AlignmentReader, fh *os.File, isBAM bool, err error) {
	fh = os.Stdin
	if file != "-" {
		fh, err = os.Open(file)
		if err != nil {
			return nil, nil, false, err
		}
	}

	br := bufio.NewReaderSize(fh, IoBuf)
	magic, err := br.Peek(2)
	if err != nil {
		if fh != os.Stdin {
			fh.Close()
		}
		return nil, nil, false, err
	}

	isBAM = magic[0] == 0x1f && magic[1] == 0x8b
	if isBAM {
		br2, err := bam.NewReader(br, 1)
		if err != nil {
			if fh != os.Stdin {
				fh.Close()
			}
			return nil, nil, false, err
		}
		return br2, fh, true, nil
	}

	sr, err := sam.NewReader(br)
	if err != nil {
		if fh != os.Stdin {
			fh.Close()
		}
		return nil, nil, false, err
	}
	return &SamReader{Reader: sr, closer: fh}, fh, false, nil
}

// DetectBarcodeFormat scans up to the first 100 records of file (SAM or
// BAM, "-" for stdin) and reports which linked-read barcode format they
// use: "haplotagging", "stlfr", "nucleotide-<LEN>", or returns an error when the
// format can't be determined -- no recognizable barcode on any sampled
// record, the file can't be opened/read, or more than one format was
// identified (input is expected to be consistently one format,
// mirroring fastq.CoreFqParser's own rule).
//
// A SAM/BAM barcode is always stored as a BX:Z tag regardless of its
// original source technology, so TELLseq and 10X -- both plain
// nucleotide strings once standardized -- are otherwise
// indistinguishable by shape alone. Rather than resolving that via the
// (expensive to load) 10X whitelist, nucleotide-shaped barcodes are
// reported as "nucleotide-<LEN>" (e.g. "nucleotide-16" for 10X, "nucleotide-18"
// for TELLseq) and left for the caller to interpret.
func DetectBarcodeFormat(file string) (string, error) {
	r, fh, isBAM, err := openAlignmentReader(file)
	if err != nil {
		return "", err
	}
	defer r.Close()
	if isBAM && fh != os.Stdin {
		defer fh.Close()
	}

	var haplotagging, stlfr, generic int
	var nucLen int

	for range 100 {
		rec, err := r.Read()
		if err != nil {
			break // EOF or a read error -- classify from whatever was sampled
		}
		bc, ok := FindBarcode(rec)
		if !ok || bc == "" {
			continue
		}
		switch {
		case isHaplotaggingBarcode(bc):
			haplotagging++
		case isStlfrBody(bc):
			stlfr++
		default:
			if l := nucleotideCoreLen(bc); l > 0 {
				generic++
				if l > nucLen {
					nucLen = l
				}
			}
		}
	}

	formatsFound := 0
	if haplotagging > 0 {
		formatsFound++
	}
	if stlfr > 0 {
		formatsFound++
	}
	if generic > 0 {
		formatsFound++
	}

	switch {
	case formatsFound > 1:
		return "", fmt.Errorf("multiple barcode formats detected, cannot proceed (multiple barcode formats may be permissble with more widespread adoption of Standard format).")
	case haplotagging > 0:
		return "haplotagging", nil
	case stlfr > 0:
		return "stlfr", nil
	case generic > 0:
		return fmt.Sprintf("nucleotide-%d", nucLen), nil
	}
	return "", fmt.Errorf("unable to detect barcode format from first 100 records")
}
