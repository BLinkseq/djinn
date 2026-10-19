// Package stats reports linked-read metrics per molecule from a
// coordinate-sorted SAM/BAM file.
//
// # Plain-language picture of what happens here
//
// Linked-read sequencing tags every DNA fragment ("molecule") with a barcode, so all
// reads that came from one molecule share a barcode. This package walks through an
// alignment file once, from the start of each chromosome to its end, and for every
// barcode it keeps a running tally of its reads: how many, where the first one
// starts, where the furthest one ends, how many bases lined up, and so on. When a
// chromosome is finished it prints one line of the table per barcode and forgets them.
//
// Reads sharing a barcode are separated into molecules where the gap between consecutive reads is
// bigger than a chosen distance (see package deconvolve, which does the separating); a molecule is
// printed as soon as the scan has moved too far past it for another read to join it.
package stats

import (
	"bufio"
	"errors"
	"io"
	"os"
	"strconv"

	"djinn/deconvolve"
	"djinn/xam"

	"github.com/biogo/hts/sam"
)

// The column names of the output table (the first line of the output).
const header = "contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n"

// ErrNotSorted is returned when the input is not coordinate sorted.
// "Coordinate sorted" means the reads are listed in order of where they sit on each
// chromosome. The whole method relies on that order, so we refuse to guess without it.
var ErrNotSorted = deconvolve.ErrNotSorted

// StatsXam writes one line per molecule of infile to stdout.
//
// Records are grouped by barcode, found with xam.FindBarcode: the BX tag, or a TELLseq/stLFR
// style barcode in the read name. A barcode is valid when the VX tag says so, or, with no VX tag,
// when the barcode itself is valid. Reads without a valid barcode are pooled into one "invalid"
// line per contig. Reads sharing a barcode are split into separate molecules (BARCODE-1, ...) when
// the gap between consecutive reads exceeds distance; distance 0 never splits.
//
// The input must be coordinate sorted. Molecules written before a sorting error are still output.
func StatsXam(infile string, distance, threads int) error {
	infile = xam.FileOrStdin(infile) // "" means "read from a pipe", so this decides file vs. stdin

	// A background worker reads the file and hands us records through a channel (a conveyor belt).
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, threads)

	// Output goes through a buffer: lines pile up in memory and are written to the screen/file in big
	// batches, which is much faster than one tiny write per line. If a write fails, the buffer remembers
	// the error and reports it when we call Flush below, so the individual writes need not be checked.
	out := bufio.NewWriterSize(os.Stdout, xam.IoBuf)
	_, _ = out.WriteString(header) // (see above: errors are reported by Flush)

	runErr := run(recChan, out, distance)
	if runErr != nil {
		// We stopped early, so the reader may still be holding unread records. Let it finish up in the
		// background instead of leaving it stuck waiting for us to take the next one.
		go func() {
			for range recChan {
			}
		}()
	}
	if err := out.Flush(); err != nil { // push out whatever is left in the buffer
		return errors.Join(runErr, err)
	}
	return runErr
}

// run reads every record from recs and writes the molecule table to w.
func run(recs <-chan *sam.Record, w io.Writer, distance int) error {
	t := newTracker(w, distance)
	for rec := range recs {
		if err := t.process(rec); err != nil {
			return err
		}
	}
	t.split.Finish() // end of the file: prints whatever molecules are still open
	return nil
}

// ── what a molecule looks like ────────────────────────────────────────────────

// molecule is the running tally for one barcode on one chromosome: a handful of counters that
// are updated as each of its reads is seen, so we never need to keep the reads themselves.
type molecule struct {
	chrom, name string // the chromosome and the barcode this tally is printed under
	start, end  int    // where the molecule begins, and the furthest point any of its reads reaches
	bp, insert  int    // total bases that lined up, and total length of the DNA fragments the reads came from
	count       int    // number of reads (a read pair counts once)
	started     bool   // false until the first read arrives
}

// add takes one ordinary (primary) read into the tally.
func (m *molecule) add(r *sam.Record) {
	end := xam.AlignedEnd(r)
	if !m.started {
		m.start, m.end, m.started = r.Pos, end, true // the first read fixes where the molecule begins
	} else {
		m.end = max(m.end, end) // the molecule reaches as far right as its furthest read
	}
	m.chrom = r.Ref.Name()
	m.bp += alignedLength(r)
	m.insert += insertSize(r)
	m.count += readCount(r)
}

// addSupplementary handles a "supplementary" alignment: an extra piece of a read that was already
// counted (for example a read that spans a rearrangement and so lines up in two places). We credit
// its bases to the molecule but leave the read count and its extent untouched. If the barcode has no
// molecule yet there is nothing to credit and it is ignored.
func (m *molecule) addSupplementary(r *sam.Record) {
	if m.started {
		m.bp += alignedLength(r)
		m.insert += insertSize(r)
	}
}

// invalidPool collects every read that has no valid barcode. They cannot be told apart, so they are not
// treated as a molecule: each chromosome gets a single "invalid" line that just totals their reads and
// aligned bases (and shows no start, end or insert).
type invalidPool struct {
	chrom     string
	bp, reads int
}

func (p *invalidPool) add(r *sam.Record) {
	p.chrom = r.Ref.Name()
	p.bp += alignedLength(r)
	p.reads += readCount(r)
}

func (p *invalidPool) asMolecule() *molecule {
	return &molecule{chrom: p.chrom, name: "invalid", bp: p.bp, count: p.reads}
}

// readCount is how many reads this record adds to a molecule's read count. A paired read
// counts only through its first mate, so a pair is counted once instead of twice.
func readCount(r *sam.Record) int {
	if r.Flags&sam.Paired == 0 || r.Flags&sam.Read1 != 0 {
		return 1
	}
	return 0
}

// ── the tracker: everything that changes while we scan the file ───────────────

// tracker holds the state of one pass through the file and decides what to do with each read.
type tracker struct {
	w    io.Writer
	line []byte // scratch space reused for each output line, so printing does not allocate each time

	split     *deconvolve.Splitter // says which molecule of its barcode each read belongs to
	molecules map[string]*molecule // the open molecule of every barcode being followed, by barcode
	invalid   *invalidPool         // reads without a valid barcode on the current chromosome, pooled together
}

func newTracker(w io.Writer, distance int) *tracker {
	t := &tracker{w: w, molecules: map[string]*molecule{}, split: deconvolve.NewSplitter(distance)}
	// Print a molecule the moment the scan has moved too far past it for another read to join it, and
	// print whatever is still open when a chromosome ends.
	t.split.OnClose(t.closeMolecule)
	t.split.OnContigEnd(t.flushAll)
	return t
}

// process deals with one record from the file.
func (t *tracker) process(rec *sam.Record) error {
	place, err := t.split.Place(rec)
	if err != nil {
		return err
	}
	switch place.Kind {
	case deconvolve.Supplementary:
		// Credit the extra piece to its barcode's open molecule, if it has one.
		if m := t.molecules[place.Barcode]; m != nil && place.Barcode != "" {
			m.addSupplementary(rec)
		}
	case deconvolve.NoBarcode:
		if t.invalid == nil {
			t.invalid = &invalidPool{}
		}
		t.invalid.add(rec)
	case deconvolve.Primary:
		m := t.molecules[place.Barcode]
		if m == nil { // first read of this molecule: start its tally, named after the barcode (plus -1, -2, ... after the first)
			m = &molecule{name: place.Barcode}
			if place.Molecule > 0 {
				m.name += "-" + strconv.Itoa(place.Molecule)
			}
			t.molecules[place.Barcode] = m
		}
		m.add(rec)
	}
	return nil // skipped, duplicate and secondary records do not count
}

// closeMolecule prints the open molecule of a barcode and forgets it.
func (t *tracker) closeMolecule(barcode string) {
	if m := t.molecules[barcode]; m != nil {
		t.write(m)
		delete(t.molecules, barcode)
	}
}

// flushAll prints every molecule still open and the pooled invalid reads, used at the end of each chromosome.
func (t *tracker) flushAll() {
	if t.invalid != nil {
		t.write(t.invalid.asMolecule())
		t.invalid = nil
	}
	for _, m := range t.molecules {
		t.write(m)
	}
	clear(t.molecules)
}

// write prints one molecule as a line of the output table.
func (t *tracker) write(m *molecule) {
	t.line = appendLine(t.line[:0], m) // rebuild the line in the reusable scratch space
	_, _ = t.w.Write(t.line)           // the buffered writer reports any error when it is flushed
}

// appendLine adds the tab-separated output line for m to dst. Building the line with strconv's
// append functions (instead of fmt) avoids creating temporary strings for every number.
func appendLine(dst []byte, m *molecule) []byte {
	inferred := m.end - m.start // the molecule's length as seen from where its reads landed

	dst = append(dst, m.chrom...)
	dst = append(dst, '\t')
	dst = append(dst, m.name...)
	for _, n := range [...]int{max(1, m.count), m.start, m.end, inferred, m.bp, m.insert} {
		dst = append(dst, '\t')
		dst = strconv.AppendInt(dst, int64(n), 10)
	}

	// Coverage: the fraction of the inferred length actually backed by aligned bases (or by fragment
	// length), capped at 1. An empty molecule (length 0) has no meaningful fraction, so it is 0.
	for _, n := range [...]int{m.bp, m.insert} {
		dst = append(dst, '\t')
		if inferred == 0 {
			dst = append(dst, "0.00000"...)
		} else {
			dst = appendCoverage(dst, n, inferred)
		}
	}
	return append(dst, '\n')
}
