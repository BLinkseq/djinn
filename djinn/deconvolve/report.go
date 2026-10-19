package deconvolve

import (
	"fmt"
	"io"
	"slices"
	"strconv"
	"strings"
)

// report is the end-of-run summary: how many molecules there were before and after splitting, and how
// many barcodes were split how many times. A "molecule" before splitting is one barcode on one chromosome
// (the same barcode on two chromosomes counts as two, since the numbering starts over on each).
type report struct {
	original  int         // barcodes (per chromosome), i.e. molecules before splitting
	splits    int         // total number of splits made
	histogram map[int]int // number of splits -> how many barcodes were split that many times
}

// add records one barcode (on one chromosome) that was split `splits` times.
func (r *report) add(splits int) {
	if r.histogram == nil {
		r.histogram = map[int]int{}
	}
	r.original++
	r.splits += splits
	r.histogram[splits]++
}

// write prints the summary as a tab-separated table:
//
//	Original_Molecules   molecules before splitting
//	New_Molecules        molecules after splitting (each split adds one)
//	Deconvolutions       Barcodes
//	0                    barcodes that were not split
//	1                    barcodes that were split once
//	...
func (r *report) write(w io.Writer) error {
	var sb strings.Builder
	row := func(name string, count int) {
		sb.WriteString(name)
		sb.WriteByte('\t')
		sb.WriteString(strconv.Itoa(count))
		sb.WriteByte('\n')
	}
	row("Original_Molecules", r.original)
	row("New_Molecules", r.original+r.splits)
	sb.WriteString("Deconvolutions\tBarcodes\n")

	counts := make([]int, 0, len(r.histogram))
	for n := range r.histogram {
		counts = append(counts, n)
	}
	slices.Sort(counts) // fewest splits first
	for _, n := range counts {
		row(strconv.Itoa(n), r.histogram[n])
	}
	_, err := fmt.Fprint(w, sb.String())
	return err
}
