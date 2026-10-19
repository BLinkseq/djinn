// Package stats reports linked-read metrics per molecule from a
// coordinate-sorted SAM/BAM file.
package stats

import (
	"bufio"
	"container/heap"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"

	"djinn/xam"

	"github.com/biogo/hts/sam"
)

const header = "contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n"

// ErrNotSorted is returned when the input is not coordinate sorted.
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
	infile = xam.FileOrStdin(infile)
	recChan, _ := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, threads)

	out := bufio.NewWriterSize(os.Stdout, xam.IoBuf)
	out.WriteString(header)
	runErr := run(recChan, out, distance)
	if err := out.Flush(); err != nil {
		return err
	}
	return runErr
}

func writeMolecule(w io.Writer, chrom, barcode string, start, end, insert, bp, count int) {
	inferred := end - start
	var covBP, covInsert float64
	if inferred != 0 {
		covBP, covInsert = coverage(bp, inferred), coverage(insert, inferred)
	}
	fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
		chrom, barcode, max(1, count), start, end, inferred, bp, insert, pyFloat(covBP), pyFloat(covInsert))
}

// ── per-barcode state ─────────────────────────────────────────────────────────

// molecule holds the running totals of the molecule being assembled.
type molecule struct {
	start, end int // end is the largest reference end seen
	lastEnd    int // reference end of the latest read, for gap detection
	bp, insert int
	count      int
}

// readCloud tracks one barcode in constant memory. The input is coordinate
// sorted, so a split between molecules is detected as reads arrive.
type readCloud struct {
	barcode, chrom string
	suffix         int
	valid          bool
	cutoff         int
	mol            *molecule
	invalidBP      int
	invalidCount   int
}

// add takes one primary read, writing a finished molecule to w when this read starts a new one.
func (c *readCloud) add(r *sam.Record, w io.Writer) {
	c.chrom = r.Ref.Name()
	bp := alignedLength(r)
	ins := insertSize(r)
	count := 0
	if r.Flags&sam.Paired == 0 || r.Flags&sam.Read1 != 0 {
		count = 1
	}

	if !c.valid {
		c.invalidBP += bp
		c.invalidCount += count
		return
	}

	start, end := r.Pos, alignedEnd(r)
	switch {
	case c.mol == nil:
		c.mol = &molecule{start, end, end, bp, ins, count}
	case c.cutoff > 0 && start-c.mol.lastEnd > c.cutoff:
		c.emit(w)
		c.suffix++
		c.mol = &molecule{start, end, end, bp, ins, count}
	default:
		c.mol.end = max(c.mol.end, end)
		c.mol.lastEnd = end
		c.mol.bp += bp
		c.mol.insert += ins
		c.mol.count += count
	}
}

// addSupplementary credits the aligned bases of a supplementary alignment to the current
// molecule without touching its count, end or lastEnd, so it can neither move the molecule
// boundary nor trigger a split. It is dropped if the barcode has no molecule yet.
func (c *readCloud) addSupplementary(r *sam.Record) {
	if c.mol != nil {
		c.mol.bp += alignedLength(r)
		c.mol.insert += insertSize(r)
	}
}

// flush writes whatever is in progress.
func (c *readCloud) flush(w io.Writer) {
	if !c.valid {
		if c.chrom != "" {
			writeMolecule(w, c.chrom, "invalid", 0, 0, 0, c.invalidBP, c.invalidCount)
		}
		return
	}
	if c.mol != nil {
		c.emit(w)
		c.mol = nil
	}
}

// lastEnd is the end of the latest read, which orders the eviction heap.
func (c *readCloud) lastEnd() int {
	if c.mol == nil {
		return 0
	}
	return c.mol.lastEnd
}

func (c *readCloud) emit(w io.Writer) {
	name := c.barcode
	if c.suffix != 0 {
		name = c.barcode + "-" + strconv.Itoa(c.suffix)
	}
	m := c.mol
	writeMolecule(w, c.chrom, name, m.start, m.end, m.insert, m.bp, m.count)
}

// ── eviction heap ─────────────────────────────────────────────────────────────

// evictEntry is a barcode and the read end it had when pushed. Entries go stale when the
// barcode gets another read; they are detected and skipped on pop.
type evictEntry struct {
	end int
	bx  string
}

type evictHeap []evictEntry

func (h evictHeap) Len() int { return len(h) }
func (h evictHeap) Less(i, j int) bool {
	return h[i].end < h[j].end || (h[i].end == h[j].end && h[i].bx < h[j].bx)
}
func (h evictHeap) Swap(i, j int) { h[i], h[j] = h[j], h[i] }
func (h *evictHeap) Push(x any)   { *h = append(*h, x.(evictEntry)) }
func (h *evictHeap) Pop() any {
	old := *h
	e := old[len(old)-1]
	*h = old[:len(old)-1]
	return e
}

// ── main loop ─────────────────────────────────────────────────────────────────

// run streams the molecule table for the records in recs to w. It returns ErrNotSorted
// for unsorted input.
func run(recs <-chan *sam.Record, w io.Writer, distance int) error {
	clouds := map[string]*readCloud{}
	// The bucket for reads without a valid barcode, kept out of clouds so no barcode can share its slot.
	var invalid *readCloud
	// Keeps a barcode's molecule numbering going if it reappears on the contig after an early eviction.
	suffixNext := map[string]int{}
	evict := &evictHeap{}
	var lastContig string
	var lastPos int

	flushAll := func() {
		if invalid != nil {
			invalid.flush(w)
			invalid = nil
		}
		for _, c := range clouds {
			c.flush(w)
		}
		clear(clouds)
		*evict = (*evict)[:0]
		clear(suffixNext)
	}

	for rec := range recs {
		if rec.Flags&(sam.Unmapped|sam.Duplicate|sam.Secondary) != 0 || len(rec.Cigar) == 0 {
			continue
		}

		chrom := rec.Ref.Name()
		if lastContig != "" && chrom != lastContig {
			flushAll()
			lastPos = 0
		}
		lastContig = chrom
		if rec.Pos < lastPos {
			return ErrNotSorted
		}
		lastPos = rec.Pos

		bx, valid := xam.FindBarcode(rec)

		if rec.Flags&sam.Supplementary != 0 {
			if c := clouds[bx]; c != nil && bx != "" {
				c.addSupplementary(rec)
			}
			continue
		}

		// A barcode whose latest read ends more than `distance` bp behind this one cannot gain
		// another read in its current molecule, so it is emitted and reclaimed now.
		for distance > 0 && evict.Len() > 0 && (*evict)[0].end+distance < rec.Pos {
			e := heap.Pop(evict).(evictEntry)
			c := clouds[e.bx]
			if c == nil || c.lastEnd() > e.end {
				continue // stale entry: the barcode was updated after this was pushed
			}
			c.flush(w)
			suffixNext[e.bx] = c.suffix + 1
			delete(clouds, e.bx)
		}

		if bx == "" || !valid {
			if invalid == nil {
				invalid = &readCloud{barcode: "invalid", chrom: chrom}
			}
			invalid.add(rec, w)
			continue
		}

		c := clouds[bx]
		if c == nil {
			c = &readCloud{barcode: bx, chrom: chrom, valid: true, suffix: suffixNext[bx], cutoff: distance}
			clouds[bx] = c
		}
		c.add(rec, w)
		if distance > 0 {
			heap.Push(evict, evictEntry{c.lastEnd(), bx})
		}
	}
	flushAll()
	return nil
}
