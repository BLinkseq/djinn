// Package deconvolve separates reads that share a barcode but came from different molecules.
//
// # Plain-language picture of what happens here
//
// A linked-read barcode is meant to say "these reads came from the same long DNA molecule".
// But a library holds many molecules, so by chance two different molecules can carry the same
// barcode. They can be told apart by where they sit on the chromosome: reads from one molecule
// huddle together, and a long empty stretch before the next read with that barcode suggests it
// came from a different molecule.
//
// This package walks through an alignment file from the start of each chromosome to its end. For
// every barcode it remembers where its latest read ended. When the next read with that barcode starts
// more than a chosen distance further along, it is taken to be a new molecule, and from then on the
// barcode is written with a number on the end: ACGT, then ACGT-1, then ACGT-2, and so on. Everything
// else about every record is left exactly as it was.
package deconvolve

import (
	"errors"
	"os"
	"strconv"

	"djinn/xam"

	"github.com/biogo/hts/sam"
)

// ErrNotSorted is returned when the input is not coordinate sorted.
// "Coordinate sorted" means the reads are listed in order of where they sit on each
// chromosome. The whole method relies on that order, so we refuse to guess without it.
var ErrNotSorted = errors.New("input alignment file does not appear to be coordinate-sorted, " +
	"which will result in incorrect output. Please sort the alignments by coordinates (e.g. samtools sort) and try again")

// DeconvolveXam reads a coordinate-sorted SAM/BAM file (or stdin for "-") and writes every record to
// stdout, with the barcode (BX tag) suffixed by "-1", "-2", ... wherever a barcode's reads are more
// than distance bases apart. The first molecule of a barcode keeps the barcode as it is. Barcode
// numbering starts over on each chromosome. Records without a valid barcode, and unmapped records,
// are written unchanged.
//
// A summary goes to stderr when the file is done: the number of molecules before and after, and
// how many barcodes were split how many times. A molecule here is one barcode on one chromosome.
//
// asSam writes SAM text instead of BAM. distance must be greater than 0.
func DeconvolveXam(infile string, distance int, asSam bool, threads int) error {
	if distance <= 0 {
		return errors.New("the distance threshold must be greater than 0")
	}
	readThread, writeThread := xam.IoThreads(threads)

	// A background worker reads the file and hands us records through a channel (a conveyor belt).
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// Note in the file's header that djinn touched it, as the other djinn commands do.
	hdr := br.Header()
	if err := hdr.AddProgram(xam.NewPG(hdr, "djinn sam deconvolve -d "+strconv.Itoa(distance)+" "+infile)); err != nil {
		return err
	}

	// Another background worker takes finished records from a second channel and writes them out.
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	t := newTracker(distance)
	var runErr error
	for rec := range recChan {
		if runErr = t.process(rec); runErr != nil {
			break
		}
		writeChan <- rec // every record is written, whether or not its barcode changed
	}
	close(writeChan)
	<-writeDone // wait until everything has really been written out

	if runErr != nil {
		// We stopped early, so the reader may still be holding unread records. Let it finish up in the
		// background instead of leaving it stuck waiting for us to take the next one.
		go func() {
			for range recChan {
			}
		}()
		return runErr
	}
	t.finishContig() // the last chromosome is still open
	return t.report.write(os.Stderr)
}

// ── what we remember about one barcode ────────────────────────────────────────

// barcodeState is everything we need to know about a barcode on the current chromosome.
type barcodeState struct {
	suffix  int // which molecule of this barcode we are in: 0 for the first, then 1, 2, ... (so also how many times it has been split)
	lastEnd int // where the barcode's most recent read ended
}

// ── the tracker: everything that changes while we scan the file ───────────────

type tracker struct {
	distance int                      // a gap bigger than this starts a new molecule
	states   map[string]*barcodeState // every barcode seen so far on the current chromosome
	report   report

	lastContig string // the chromosome of the previous read
	lastPos    int    // where the previous read started, to check the file really is sorted

	// Where the most recent ordinary read started. It tells us which barcodes have gone quiet; see label.
	lastPrimaryPos int
}

func newTracker(distance int) *tracker {
	return &tracker{distance: distance, states: map[string]*barcodeState{}}
}

// process looks at one record and, if its barcode belongs to a later molecule, changes the record's BX tag.
func (t *tracker) process(rec *sam.Record) error {
	// Reads that are not aligned anywhere carry no position, so there is nothing to place them by.
	if rec.Flags&sam.Unmapped != 0 || len(rec.Cigar) == 0 {
		return nil
	}
	chrom := rec.Ref.Name()

	// Duplicates and secondary alignments are not used to decide where molecules begin and end (so they
	// are also not checked for sort order), but they still get the label of the molecule they sit in.
	if rec.Flags&(sam.Duplicate|sam.Secondary) != 0 {
		if chrom == t.lastContig { // before the first usable read of a new chromosome there is nothing to compare with
			t.relabel(rec)
		}
		return nil
	}

	// Moving to a new chromosome ends everything in progress: the numbering starts over.
	if t.lastContig != "" && chrom != t.lastContig {
		t.finishContig()
		t.lastPos, t.lastPrimaryPos = 0, 0
	}
	t.lastContig = chrom

	// The file must be in left-to-right order within each chromosome.
	if rec.Pos < t.lastPos {
		return ErrNotSorted
	}
	t.lastPos = rec.Pos

	// An extra piece of a read that was already counted ("supplementary") never starts a molecule;
	// it just gets the label of the molecule its barcode is in.
	if rec.Flags&sam.Supplementary != 0 {
		t.relabel(rec)
		return nil
	}

	t.lastPrimaryPos = rec.Pos
	bx, valid := xam.FindBarcode(rec)
	if bx == "" || !valid {
		return nil // no usable barcode, nothing to split
	}

	t.advance(bx, rec)
	setSuffixed(rec, bx, t.states[bx].suffix)
	return nil
}

// advance takes an ordinary read with a valid barcode into account: it decides whether this read starts
// a new molecule of its barcode (and if so counts the split), and remembers where the read ended.
func (t *tracker) advance(bx string, rec *sam.Record) {
	st := t.states[bx]
	if st == nil { // first read of this barcode on this chromosome: its first molecule
		t.states[bx] = &barcodeState{lastEnd: xam.AlignedEnd(rec)}
		return
	}
	// The gap is measured from where the previous read ended (not from the furthest point the molecule
	// reached), so one very long read does not hide a gap after the shorter reads that follow it.
	if rec.Pos-st.lastEnd > t.distance {
		st.suffix++
	}
	st.lastEnd = xam.AlignedEnd(rec)
}

// relabel gives a supplementary, duplicate or secondary record the label of its barcode's molecule.
func (t *tracker) relabel(rec *sam.Record) {
	bx, _ := xam.FindBarcode(rec)
	if bx == "" {
		return
	}
	setSuffixed(rec, bx, t.label(bx))
}

// label is the molecule number a record that does not itself move the barcode along should carry.
//
// If the barcode has not been seen yet, it is the first molecule (0). If it has, the record belongs to the
// barcode's current molecule, unless that molecule has already been closed: once the scan has moved on
// (by an ordinary read) more than the distance past the barcode's last read, that molecule can take no
// more reads, so the record belongs to the molecule that would come next, number+1. The numbering is
// what lets a later step (such as `djinn sam stats`) count a record towards a molecule only while that
// molecule is still open.
func (t *tracker) label(bx string) int {
	st := t.states[bx]
	if st == nil {
		return 0
	}
	if t.lastPrimaryPos > st.lastEnd+t.distance {
		return st.suffix + 1
	}
	return st.suffix
}

// finishContig adds the chromosome's barcodes to the summary and forgets them, ready for the next chromosome.
func (t *tracker) finishContig() {
	for _, st := range t.states {
		t.report.add(st.suffix)
	}
	clear(t.states)
}

// setSuffixed writes the barcode with "-n" on the end into the record's BX tag. Number 0 is the barcode's
// first molecule, which keeps the barcode exactly as it was (the record is left untouched).
func setSuffixed(rec *sam.Record, bx string, n int) {
	if n == 0 {
		return
	}
	suffixed := bx + "-" + strconv.Itoa(n)
	xam.SetBxString(rec, &suffixed) // replaces the BX tag, or adds one if the barcode came from the read name
}
