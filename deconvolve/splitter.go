package deconvolve

import (
	"container/heap"

	"github.com/blinkseq/djinn/xam"

	"github.com/biogo/hts/sam"
)

// Kind says what a record is, as far as telling molecules apart goes.
type Kind int

const (
	// Skipped records are not aligned anywhere (or have no alignment recipe), so they cannot be placed.
	Skipped Kind = iota
	// Bystander records are duplicates or secondary alignments. They never start or extend a molecule,
	// but they have a molecule number (if their barcode is being followed).
	Bystander
	// Supplementary records are extra pieces of a read that was already counted. They never start or
	// extend a molecule either, but have a molecule number.
	Supplementary
	// NoBarcode records are ordinary reads without a valid barcode: nothing to tell apart.
	NoBarcode
	// Primary records are ordinary reads with a valid barcode. They are what molecules are made of.
	Primary
)

// Placement is the answer for one record.
type Placement struct {
	Kind     Kind
	Barcode  string // the record's barcode; "" if it has none (or none that is being followed)
	Molecule int    // which molecule of that barcode the record belongs to: 0 for the first, then 1, 2, ...
}

// barcodeState is everything we need to know about a barcode on the current chromosome.
type barcodeState struct {
	suffix  int // which molecule of this barcode we are in (so also how many times it has been split)
	lastEnd int // where the barcode's most recent read ended
}

// Splitter follows a coordinate-sorted stream of records and says which molecule of its barcode each
// one belongs to. Two reads with the same barcode belong to different molecules when the gap between the
// end of the previous read and the start of the next is bigger than the distance.
//
// The gap is measured from where the previous read ended (not from the furthest point the molecule has
// reached), so one very long read does not hide a gap after the shorter reads that follow it.
//
// By default every barcode on a chromosome is remembered until the chromosome ends. A caller that wants to
// handle each molecule as soon as it is complete (and free it) registers OnClose, which also turns on the
// deadline list (see deadlines.go) that finds molecules that can take no more reads.
type Splitter struct {
	distance int // a gap bigger than this starts a new molecule; 0 means never split

	onClose     func(barcode string) // called when a barcode's molecule can take no more reads, if set
	onContigEnd func()               // called when a chromosome ends, if set

	states map[string]*barcodeState // every barcode being followed on the current chromosome

	// If a barcode is closed early and then shows up again further along the same chromosome, it
	// must carry on numbering rather than restart. This remembers the next number for each such barcode.
	suffixNext map[string]int

	deadlines evictHeap // barcodes by how soon they are overdue; used only when onClose is set

	lastContig string // the chromosome of the previous read
	lastPos    int    // where the previous read started, to check the file really is sorted

	// Where the most recent ordinary read started. It tells us which barcodes have gone quiet: a
	// molecule whose last read ended more than `distance` before this is closed. See label.
	lastPrimaryPos int
}

// NewSplitter returns a Splitter for the given distance (0 means never split).
func NewSplitter(distance int) *Splitter {
	return &Splitter{distance: distance, states: map[string]*barcodeState{}, suffixNext: map[string]int{}}
}

// OnClose registers f to be called for each molecule (named by its barcode) as soon as it can take no more
// reads: when the scan has moved more than the distance past its last read. Molecules still open when a chromosome
// ends are not reported here (see OnContigEnd). With a distance of 0 nothing is ever closed early.
func (s *Splitter) OnClose(f func(barcode string)) { s.onClose = f }

// OnContigEnd registers f to be called when a chromosome has been completely read (and at Finish), before the
// Splitter forgets that chromosome's barcodes, so f may still call Barcodes.
func (s *Splitter) OnContigEnd(f func()) { s.onContigEnd = f }

// Barcodes calls f for every barcode being followed right now, with the number of times it has been split
// (so molecules-1). Barcodes closed early are no longer included.
func (s *Splitter) Barcodes(f func(barcode string, splits int)) {
	for bx, st := range s.states {
		f(bx, st.suffix)
	}
}

// Finish must be called after the last record: it ends the last chromosome.
func (s *Splitter) Finish() { s.endContig() }

// Place says what the record is and which molecule it belongs to. It returns ErrNotSorted
// if the records are not in coordinate order.
func (s *Splitter) Place(rec *sam.Record) (Placement, error) {
	// Reads that are not aligned anywhere carry no position, so there is nothing to place them by.
	if rec.Flags&sam.Unmapped != 0 || len(rec.Cigar) == 0 {
		return Placement{Kind: Skipped}, nil
	}
	chrom := rec.Ref.Name()

	// Duplicates and secondary alignments are not used to decide where molecules begin and end (so they
	// are also not checked for sort order), but they still have the molecule number of their barcode.
	if rec.Flags&(sam.Duplicate|sam.Secondary) != 0 {
		p := Placement{Kind: Bystander}
		if chrom == s.lastContig { // before the first usable read of a new chromosome there is nothing to compare with
			if bx, _ := xam.FindBarcode(rec); bx != "" {
				p.Barcode, p.Molecule = bx, s.label(bx)
			}
		}
		return p, nil
	}

	// Moving to a new chromosome ends everything in progress: the numbering starts over.
	if s.lastContig != "" && chrom != s.lastContig {
		s.endContig()
		s.lastPos, s.lastPrimaryPos = 0, 0
	}
	s.lastContig = chrom

	// The file must be in left-to-right order within each chromosome.
	if rec.Pos < s.lastPos {
		return Placement{}, ErrNotSorted
	}
	s.lastPos = rec.Pos

	bx, valid := xam.FindBarcode(rec)

	// An extra piece of a read that was already counted never starts a molecule; it just has the number
	// of the molecule its barcode is in.
	if rec.Flags&sam.Supplementary != 0 {
		return Placement{Kind: Supplementary, Barcode: bx, Molecule: s.label(bx)}, nil
	}

	s.lastPrimaryPos = rec.Pos
	s.closeOverdue(rec.Pos) // first, report any molecule the file has moved too far past
	if bx == "" || !valid {
		return Placement{Kind: NoBarcode}, nil
	}
	return Placement{Kind: Primary, Barcode: bx, Molecule: s.advance(bx, rec)}, nil
}

// advance takes an ordinary read with a valid barcode into account and returns its molecule number: it decides
// whether this read starts a new molecule of its barcode, and remembers where the read ended.
func (s *Splitter) advance(bx string, rec *sam.Record) int {
	st := s.states[bx]
	switch {
	case st == nil: // first read of this barcode (on this chromosome, or since it was closed early)
		st = &barcodeState{suffix: s.suffixNext[bx]}
		s.states[bx] = st
		delete(s.suffixNext, bx)
	case s.distance > 0 && rec.Pos-st.lastEnd > s.distance:
		// The gap since the previous read is too big: this read belongs to a different molecule that
		// happens to carry the same barcode. (With OnClose set, such a barcode was already closed by
		// closeOverdue, so this is only reached when molecules are kept until the chromosome ends.)
		st.suffix++
	}
	st.lastEnd = xam.AlignedEnd(rec)

	if s.onClose != nil && s.distance > 0 {
		// The barcode's latest read now ends somewhere new, so add a fresh line to the deadline list.
		// (The barcode's older line is left behind as a stale entry; see deadlines.go.)
		heap.Push(&s.deadlines, evictEntry{end: st.lastEnd, bx: bx})
	}
	return st.suffix
}

// label is the molecule number a record that does not itself move the barcode along (supplementary,
// duplicate, secondary) should carry.
//
// If the barcode is not being followed, it is the number its next molecule would have (0 if it has none yet).
// If it is, the record belongs to its current molecule, unless that molecule has already been closed: once the
// scan has moved on (by an ordinary read) more than the distance past the barcode's last read, that molecule can
// take no more reads, so the record belongs to the molecule that would come next, number+1. The numbering is
// what lets a later step count a record towards a molecule only while that molecule is still open.
func (s *Splitter) label(bx string) int {
	st := s.states[bx]
	if st == nil {
		return s.suffixNext[bx]
	}
	if s.distance > 0 && s.lastPrimaryPos > st.lastEnd+s.distance {
		return st.suffix + 1
	}
	return st.suffix
}

// closeOverdue reports and forgets every barcode whose latest read ended more than `distance` bases before
// `pos`, the place the scan has now reached. Such a barcode cannot gain another read in its current molecule,
// so keeping it in memory any longer would be wasted effort.
func (s *Splitter) closeOverdue(pos int) {
	if s.onClose == nil {
		return
	}
	// Only look at the front of the deadline list: it holds the barcode closest to being overdue.
	for s.distance > 0 && len(s.deadlines) > 0 && s.deadlines[0].end+s.distance < pos {
		entry := heap.Pop(&s.deadlines).(evictEntry) // take the front entry off the list
		st := s.states[entry.bx]

		// Is this entry stale? The barcode may be gone already, or it may have received a newer read
		// since this entry was written (so a newer, further-right entry for it is still on the list).
		// Either way this old entry is out of date, so discard it and move on.
		if st == nil || st.lastEnd > entry.end {
			continue
		}

		// The entry is current and overdue: report the molecule and stop following the barcode.
		s.onClose(entry.bx)
		s.suffixNext[entry.bx] = st.suffix + 1 // if the barcode returns later, number the next molecule onward
		delete(s.states, entry.bx)
	}
}

// endContig tells the caller the chromosome is done and forgets its barcodes, ready for the next chromosome.
func (s *Splitter) endContig() {
	if s.onContigEnd != nil {
		s.onContigEnd()
	}
	clear(s.states)
	clear(s.suffixNext)
	s.deadlines = s.deadlines[:0] // the deadline list is for one chromosome only
}
