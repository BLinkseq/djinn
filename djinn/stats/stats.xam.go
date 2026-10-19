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
// molecule is complete it prints one line of the table and forgets it.
package stats

import (
	"bufio"
	"container/heap"
	"errors"
	"io"
	"os"
	"strconv"

	"djinn/xam"

	"github.com/biogo/hts/sam"
)

// The column names of the output table (the first line of the output).
const header = "contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n"

// ErrNotSorted is returned when the input is not coordinate sorted.
// "Coordinate sorted" means the reads are listed in order of where they sit on each
// chromosome. The whole method relies on that order, so we refuse to guess without it.
var ErrNotSorted = errors.New("input alignment file does not appear to be coordinate-sorted, " +
	"which will result in incorrect output statistics. Please sort the alignments by coordinates (e.g. samtools sort) and try again")

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
	t.flushAll() // end of the file: print whatever molecules are still open
	return nil
}

// ── what a finished molecule looks like ───────────────────────────────────────

// molecule is the running tally for one molecule: a handful of counters that are
// updated as each of its reads is seen, so we never need to keep the reads themselves.
type molecule struct {
	start, end int // where the molecule begins, and the furthest point any of its reads reaches
	lastEnd    int // where the most recent read ended (used to measure the gap to the next read)
	bp, insert int // total bases that lined up, and total length of the DNA fragments the reads came from
	count      int // number of reads (a read pair counts once)
}

// finishedMolecule is a completed molecule together with the label it is printed under.
type finishedMolecule struct {
	chrom, name string // the chromosome, and the barcode (with -1, -2, ... for later molecules)
	molecule
}

// ── one barcode's state ───────────────────────────────────────────────────────

// readCloud follows a single barcode. Because the file is sorted, a barcode's reads arrive in
// left-to-right order, so we can tell when one molecule ends and the next begins just by looking
// at the gap between consecutive reads. Only one molecule per barcode is ever held at a time.
type readCloud struct {
	barcode string
	chrom   string
	suffix  int      // 0 for the first molecule of this barcode, then 1, 2, ... for later ones
	cutoff  int      // a gap bigger than this starts a new molecule (0 = never split)
	mol     molecule // the molecule being built right now
	open    bool     // false until the first read arrives, and again after the molecule is printed
}

// add takes one ordinary (primary) read. If this read starts a new molecule, the previous
// one is complete and is returned (second result true) so the caller can print it.
func (c *readCloud) add(r *sam.Record) (done finishedMolecule, finished bool) {
	c.chrom = r.Ref.Name()
	// A fresh molecule built from just this read. It is used if this read starts a new molecule.
	start, end := r.Pos, alignedEnd(r)
	fresh := molecule{
		start: start, end: end, lastEnd: end,
		bp: alignedLength(r), insert: insertSize(r), count: readCount(r),
	}

	switch {
	case !c.open:
		// First read of this barcode on this stretch: start a molecule.
		c.mol, c.open = fresh, true

	case c.cutoff > 0 && start-c.mol.lastEnd > c.cutoff:
		// The gap since the previous read is too big: this read belongs to a different molecule that
		// happens to carry the same barcode. Hand back the old molecule and start a new one.
		done, finished = c.finish(), true
		c.suffix++
		c.mol = fresh

	default:
		// Same molecule: fold this read into the running tally.
		c.mol.end = max(c.mol.end, end) // the molecule reaches as far right as its furthest read
		c.mol.lastEnd = end             // remember where this read ended, for measuring the next gap
		c.mol.bp += fresh.bp
		c.mol.insert += fresh.insert
		c.mol.count += fresh.count
	}
	return done, finished
}

// addSupplementary handles a "supplementary" alignment: an extra piece of a read that was already
// counted (for example a read that spans a rearrangement and so lines up in two places). We credit
// its bases to the open molecule but leave the read count and the gap measurement untouched, so
// it can neither inflate the read count nor cause a split. If the barcode has no open molecule
// yet, there is nothing to credit and it is ignored.
func (c *readCloud) addSupplementary(r *sam.Record) {
	if c.open {
		c.mol.bp += alignedLength(r)
		c.mol.insert += insertSize(r)
	}
}

// flush closes the open molecule (if any) and returns it for printing.
func (c *readCloud) flush() (done finishedMolecule, finished bool) {
	if !c.open {
		return done, false
	}
	c.open = false
	return c.finish(), true
}

// finish labels the molecule being built: the barcode, plus "-1", "-2", ... after the first one.
func (c *readCloud) finish() finishedMolecule {
	name := c.barcode
	if c.suffix != 0 {
		name += "-" + strconv.Itoa(c.suffix)
	}
	return finishedMolecule{chrom: c.chrom, name: name, molecule: c.mol}
}

// lastEnd is where this barcode's most recent read ended (0 if nothing is open). It is the key
// the "deadline list" below is ordered by.
func (c *readCloud) lastEnd() int {
	if !c.open {
		return 0
	}
	return c.mol.lastEnd
}

// ── reads without a usable barcode ────────────────────────────────────────────

// invalidPool collects every read that has no valid barcode. They cannot be told apart, so
// they are not split into molecules: each chromosome gets a single "invalid" line that just
// totals their reads and aligned bases.
type invalidPool struct {
	chrom     string
	bp, reads int
}

func (p *invalidPool) add(r *sam.Record) {
	p.chrom = r.Ref.Name()
	p.bp += alignedLength(r)
	p.reads += readCount(r)
}

func (p *invalidPool) finished() finishedMolecule {
	return finishedMolecule{chrom: p.chrom, name: "invalid", molecule: molecule{bp: p.bp, count: p.reads}}
}

// readCount is how many reads this record adds to a molecule's read count. A paired read
// counts only through its first mate, so a pair is counted once instead of twice.
func readCount(r *sam.Record) int {
	if r.Flags&sam.Paired == 0 || r.Flags&sam.Read1 != 0 {
		return 1
	}
	return 0
}

// ── the deadline list (a heap) ────────────────────────────────────────────────
//
// The problem: a barcode's molecule is only finished once we know no more of its reads are
// coming, and we can only know that by noticing that the file has moved far past it. Waiting
// until the end of the chromosome would mean holding every barcode in memory the whole time.
//
// The idea: keep a to-do list of barcodes, ordered by "whose most recent read ended furthest
// to the left". That barcode is the one closest to being finished. Each time we move to a new
// read, we only need to look at the top of the list: if the file has moved more than
// `distance` bases past the top barcode's last read, its molecule is complete (any later read
// would have been too far away to join it anyway), so we print it, drop it from memory, and
// check the next barcode on the list. If the top barcode is not overdue, none of the others
// are either, because the list is ordered, so we stop looking immediately.
//
// A "heap" is the standard tool for such a list: it always hands back the smallest item
// first, and adding or removing an item is quick (it does not re-sort everything). Go's
// container/heap does the work, as long as we teach it five small things about our list:
// how long it is, which of two entries comes first, how to swap two entries, and how to add
// and remove the last one.
//
// One subtlety: when a barcode gets another read, its position in the list should change, but
// a heap cannot cheaply edit an item in the middle. So we do not edit: we add a new entry
// with the barcode's new position and leave the old entry where it is. The old entry is now
// "stale" (out of date). When a stale entry reaches the top, we notice it no longer matches
// the barcode's real position and simply throw it away. This trick is called lazy deletion.

// evictEntry is one line on the deadline list: a barcode and where its latest read ended
// at the moment the line was written.
type evictEntry struct {
	end int
	bx  string
}

// evictHeap is the deadline list. It lists the entries in an order where the entry with the
// smallest `end` (the one furthest left) is always first.
type evictHeap []evictEntry

func (h evictHeap) Len() int { return len(h) } // how many entries there are

// Less says which of two entries goes first: the one whose read ended further left. If they
// tie, the barcode name breaks the tie, so the order is always predictable.
func (h evictHeap) Less(i, j int) bool {
	return h[i].end < h[j].end || (h[i].end == h[j].end && h[i].bx < h[j].bx)
}

func (h evictHeap) Swap(i, j int) { h[i], h[j] = h[j], h[i] } // exchange two entries' places

// Push and Pop are called by container/heap, not by us directly: heap.Push(&h, x) puts x at
// the end and then shuffles it up into its right place, and heap.Pop(&h) moves the first entry
// to the end, re-orders the rest, and then calls Pop to cut the last one off and return it.
func (h *evictHeap) Push(x any) { *h = append(*h, x.(evictEntry)) }
func (h *evictHeap) Pop() any {
	old := *h
	last := old[len(old)-1]
	*h = old[:len(old)-1]
	return last
}

// ── the tracker: everything that changes while we scan the file ───────────────

// tracker holds the state of one pass through the file and decides what to do with each read.
type tracker struct {
	w        io.Writer
	line     []byte // scratch space reused for each output line, so printing does not allocate each time
	distance int    // the gap (in bases) that splits a barcode's reads into separate molecules; 0 = never

	clouds  map[string]*readCloud // every barcode currently being followed, by barcode
	invalid *invalidPool          // reads without a valid barcode, on the current chromosome

	// If a barcode is closed early and then shows up again further along the same chromosome, it
	// must carry on numbering (BARCODE-1, -2, ...) rather than restart. This remembers the next
	// number for each such barcode.
	suffixNext map[string]int

	deadlines evictHeap // the deadline list described above

	lastContig string // the chromosome of the previous read
	lastPos    int    // where the previous read started, to check the file really is sorted
}

func newTracker(w io.Writer, distance int) *tracker {
	return &tracker{
		w:          w,
		distance:   distance,
		clouds:     map[string]*readCloud{},
		suffixNext: map[string]int{},
	}
}

// process deals with one record from the file.
func (t *tracker) process(rec *sam.Record) error {
	// Skip reads that should not count: not aligned anywhere, flagged as PCR duplicates, secondary
	// (alternative) alignments of a read, or with no alignment recipe (CIGAR) to measure.
	if rec.Flags&(sam.Unmapped|sam.Duplicate|sam.Secondary) != 0 || len(rec.Cigar) == 0 {
		return nil
	}

	// Moving to a new chromosome ends everything in progress: print all open molecules and start fresh.
	chrom := rec.Ref.Name()
	if t.lastContig != "" && chrom != t.lastContig {
		t.flushAll()
		t.lastPos = 0
	}
	t.lastContig = chrom

	// The file must be in left-to-right order within each chromosome.
	if rec.Pos < t.lastPos {
		return ErrNotSorted
	}
	t.lastPos = rec.Pos

	// Work out the read's barcode and whether it is valid (see StatsXam for the rules).
	bx, valid := xam.FindBarcode(rec)

	if rec.Flags&sam.Supplementary != 0 {
		if c := t.clouds[bx]; c != nil && bx != "" {
			c.addSupplementary(rec)
		}
		return nil // supplementary pieces never start molecules, so we are done with this record
	}

	t.closeOverdue(rec.Pos) // first, print any molecule the file has moved too far past

	if bx == "" || !valid {
		if t.invalid == nil {
			t.invalid = &invalidPool{}
		}
		t.invalid.add(rec)
		return nil
	}
	t.addPrimary(rec, chrom, bx)
	return nil
}

// addPrimary adds an ordinary read with a valid barcode to that barcode's molecule.
func (t *tracker) addPrimary(rec *sam.Record, chrom, bx string) {
	c := t.clouds[bx]
	if c == nil { // first time we see this barcode here: start following it
		c = &readCloud{barcode: bx, chrom: chrom, suffix: t.suffixNext[bx], cutoff: t.distance}
		t.clouds[bx] = c
	}
	if done, finished := c.add(rec); finished {
		t.write(done)
	}
	if t.distance > 0 {
		// The barcode's latest read now ends somewhere new, so add a fresh line to the deadline list.
		// (The barcode's older line is left behind as a stale entry; see the heap notes above.)
		heap.Push(&t.deadlines, evictEntry{end: c.lastEnd(), bx: bx})
	}
}

// closeOverdue prints and forgets every barcode whose latest read ended more than `distance`
// bases before `pos`, the place the scan has now reached. Such a barcode cannot gain another
// read in its current molecule, so keeping it in memory any longer would be wasted effort.
func (t *tracker) closeOverdue(pos int) {
	// Only look at the front of the deadline list: it holds the barcode closest to being overdue.
	for t.distance > 0 && len(t.deadlines) > 0 && t.deadlines[0].end+t.distance < pos {
		entry := heap.Pop(&t.deadlines).(evictEntry) // take the front entry off the list
		c := t.clouds[entry.bx]

		// Is this entry stale? The barcode may be gone already, or it may have received a newer read
		// since this entry was written (so a newer, further-right entry for it is still on the list).
		// Either way this old entry is out of date, so discard it and move on.
		if c == nil || c.lastEnd() > entry.end {
			continue
		}

		// The entry is current and overdue: print the molecule and stop following the barcode.
		if done, finished := c.flush(); finished {
			t.write(done)
		}
		t.suffixNext[entry.bx] = c.suffix + 1 // if the barcode returns later, number the next molecule onward
		delete(t.clouds, entry.bx)
	}
}

// flushAll prints every molecule still open and resets the state, used at the end of each
// chromosome and at the end of the file.
func (t *tracker) flushAll() {
	if t.invalid != nil {
		t.write(t.invalid.finished())
		t.invalid = nil
	}
	for _, c := range t.clouds {
		if done, finished := c.flush(); finished {
			t.write(done)
		}
	}
	clear(t.clouds)
	t.deadlines = t.deadlines[:0] // the deadline list is for one chromosome only
	clear(t.suffixNext)
}

// write prints one finished molecule as a line of the output table.
func (t *tracker) write(f finishedMolecule) {
	t.line = appendLine(t.line[:0], f) // rebuild the line in the reusable scratch space
	_, _ = t.w.Write(t.line)           // the buffered writer reports any error when it is flushed
}

// appendLine adds the tab-separated output line for f to dst. Building the line with strconv's
// append functions (instead of fmt) avoids creating temporary strings for every number.
func appendLine(dst []byte, f finishedMolecule) []byte {
	inferred := f.end - f.start // the molecule's length as seen from where its reads landed

	// Coverage: the fraction of the inferred length actually backed by aligned bases (or by fragment
	// length), capped at 1. An empty molecule (length 0) has no meaningful fraction, so it is 0.
	var covBP, covInsert float64
	if inferred != 0 {
		covBP, covInsert = coverage(f.bp, inferred), coverage(f.insert, inferred)
	}

	dst = append(dst, f.chrom...)
	dst = append(dst, '\t')
	dst = append(dst, f.name...)
	for _, n := range [...]int{max(1, f.count), f.start, f.end, inferred, f.bp, f.insert} {
		dst = append(dst, '\t')
		dst = strconv.AppendInt(dst, int64(n), 10)
	}
	dst = append(dst, '\t')
	dst = appendPyFloat(dst, covBP)
	dst = append(dst, '\t')
	dst = appendPyFloat(dst, covInsert)
	return append(dst, '\n')
}
