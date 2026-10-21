package stats

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/blinkseq/djinn/fastq"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// The linked-read technologies StatsFQ understands. Each hides its barcode in a different place
// (haplotagging: a BX:Z: tag in the read's comment; stLFR and TELLseq: the end of the read name)
// and each has its own idea of what makes a barcode bad, so each gets its own scheme below.
const (
	Haplotagging = "haplotagging"
	Stlfr        = "stlfr"
	Tellseq      = "tellseq"
	TenX         = "10x"
)

// StatsFQ prints a summary of the barcodes in one FASTQ file to stdout:
//
//	Reads             total sequences in the file
//	Barcodes          reads that carry a linked-read barcode
//	Barcodes_Valid    barcodes that are valid (see fastq.CoreFq: the VX tag if there is one, otherwise the barcode itself)
//	Barcodes_Invalid  barcodes that are not valid
//
// followed by one line per barcode position that can go wrong, with the number of barcodes that
// were bad at that position: haplotagging A, C, B, D (a segment reading 00, e.g. A00), stLFR 1, 2, 3
// (a part reading 0), and TELLseq 1 to 18 (an N at that base). A barcode with several bad
// positions adds to each of them.
//
// platform is "haplotagging", "stlfr" or "tellseq" (any capitalization); "" works it out from the
// first 100 reads, the way the other djinn FASTQ commands do.
func StatsFQ(fq, platform string, bcMap string) error {
	var parse func(rec *fastx.Record) (fastq.CoreFq, bool)
	var sch scheme

	if bcMap != "" {
		var bclist map[string]struct{}
		bclist, err := fastq.Load10X(bcMap)
		if err != nil {
			return err
		}
		parse = func(rec *fastx.Record) (fastq.CoreFq, bool) { return fastq.TenXR12Corefq(rec, bclist) }
		sch = tenxScheme
	}

	parse, sch, err := chooseParser(fq, platform)
	if err != nil {
		return err
	}

	seq.ValidateSeq = false // reads can hold N and other ambiguity codes; nothing here needs them checked
	reader, err := fastx.NewReader(seq.DNA, fq, "")
	if err != nil {
		return fmt.Errorf("opening %s: %w", fq, err)
	}
	defer reader.Close()

	var reads, barcodes, valid int
	bad := make([]int, sch.slots) // how many barcodes were bad at each position

	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("reading %s: %w", fq, err)
		}
		reads++

		core, ok := parse(rec)
		if !ok {
			continue // this read has no barcode
		}
		barcodes++
		if core.Valid {
			valid++
		}
		bad = sch.mark(stripSuffix(core.Barcode), bad)
	}

	_, err = io.WriteString(os.Stdout, summary(sch, reads, barcodes, valid, bad))
	return err
}

// chooseParser returns the record parser for the requested (or detected) technology,
// together with its scheme for counting bad positions.
func chooseParser(fq, platform string) (func(*fastx.Record) (fastq.CoreFq, bool), scheme, error) {
	if platform == "" {
		parse, detected, err := fastq.CoreFqParser(fq)
		if err != nil {
			return nil, scheme{}, err
		}
		if parse == nil {
			return nil, scheme{}, fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fq)
		}
		platform = detected
	}
	switch strings.ToLower(platform) {
	case Haplotagging:
		return fastq.Haplotag2Corefq, haplotaggingScheme, nil
	case Stlfr:
		return fastq.Stlfr2Corefq, stlfrScheme, nil
	case Tellseq:
		return fastq.Tellseq2Corefq, tellseqScheme, nil
	}
	return nil, scheme{}, fmt.Errorf("unknown platform %q, expected haplotagging, stlfr or tellseq", platform)
}

// stripSuffix drops an optional "-something" tail (like a GEM well number, "-1") from a barcode,
// leaving only the part that carries the technology's own positions.
func stripSuffix(bc []byte) []byte {
	if before, _, found := bytes.Cut(bc, []byte{'-'}); found {
		return before
	}
	return bc
}

// summary formats the table. It is only a couple of dozen short lines, so it is built in memory and
// written in one go.
func summary(sch scheme, reads, barcodes, valid int, bad []int) string {
	var sb strings.Builder
	row := func(name string, count int) {
		sb.WriteString(name)
		sb.WriteByte('\t')
		sb.WriteString(strconv.Itoa(count))
		sb.WriteByte('\n')
	}
	row("Reads", reads)
	row("Barcodes", barcodes)
	row("Barcodes_Valid", valid)
	row("Barcodes_Invalid", barcodes-valid)
	for i, n := range bad {
		row(sch.label(i), n)
	}
	return sb.String()
}

// ── what counts as "bad at a position" for each technology ────────────────────

// scheme describes the positions of one technology's barcode that can be bad.
type scheme struct {
	slots int                              // number of positions, always reported even when no barcode was bad there
	label func(i int) string               // the name printed for position i (counting from 0)
	mark  func(bc []byte, bad []int) []int // adds one barcode's bad positions to bad (and may lengthen it)
}

// Haplotagging barcodes look like A01C02B03D04: four segments, each a letter and two digits. A segment
// reading 00 (A00, C00, B00 or D00) means that part of the barcode could not be read.
var haplotaggingScheme = scheme{
	slots: 4,
	label: func(i int) string { return string("ACBD"[i]) }, // the order the barcode is written and printed in
	mark: func(bc []byte, bad []int) []int {
		for i := 1; i+1 < len(bc); i++ {
			if bc[i] != '0' || bc[i+1] != '0' {
				continue
			}
			switch bc[i-1] { // the letter just before the "00" says which segment failed
			case 'A', 'a':
				bad[0]++
			case 'C', 'c':
				bad[1]++
			case 'B', 'b':
				bad[2]++
			case 'D', 'd':
				bad[3]++
			}
		}
		return bad
	},
}

// stLFR barcodes look like 12_345_6: three numbers joined by underscores. A part that is exactly 0
// is a part that could not be read.
var stlfrScheme = scheme{
	slots: 3,
	label: func(i int) string { return strconv.Itoa(i + 1) },
	mark: func(bc []byte, bad []int) []int {
		part := 0 // which of the three numbers we are in
		start := 0
		for i := 0; i <= len(bc); i++ {
			if i < len(bc) && bc[i] != '_' {
				continue
			}
			if part < len(bad) && string(bc[start:i]) == "0" { // Go does not allocate for this comparison
				bad[part]++
			}
			part++
			start = i + 1
		}
		return bad
	},
}

// TELLseq barcodes are 18 bases. An N is a base that could not be called, and we count it against
// the position (1 to 18) where it sits. Longer barcodes are not an error: the table simply gets longer.
var tellseqScheme = scheme{
	slots: 18,
	label: func(i int) string { return strconv.Itoa(i + 1) },
	mark: func(bc []byte, bad []int) []int {
		for len(bad) < len(bc) {
			bad = append(bad, 0)
		}
		for i, b := range bc {
			if b == 'N' || b == 'n' {
				bad[i]++
			}
		}
		return bad
	},
}

// 10X barcodes are 16 bases. An N is a base that could not be called, and we count it against
// the position (1 to 16) where it sits. Longer barcodes are not an error: the table simply gets longer.
var tenxScheme = scheme{
	slots: 16,
	label: func(i int) string { return strconv.Itoa(i + 1) },
	mark: func(bc []byte, bad []int) []int {
		for len(bad) < len(bc) {
			bad = append(bad, 0)
		}
		for i, b := range bc {
			if b == 'N' || b == 'n' {
				bad[i]++
			}
		}
		return bad
	},
}
