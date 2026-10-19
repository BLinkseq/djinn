package deconvolve

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
// would have been too far away to join it anyway), so we report it, drop it from memory, and
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
