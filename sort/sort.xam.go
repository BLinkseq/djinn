package sort

import (
	"bytes"
	"errors"
	"github.com/blinkseq/djinn/xam"
	"io"
	"os"
	"slices"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// SortByBX streams inPath, spills BX:Z-sorted chunks in parallel, k-way merges into outPath.
func SortByBX(infile, outPath, tmpDir string, threads int, asSam bool) error {
	// ── open reader -------------───────────────────────────────────────────────
	// Batched: the reader goroutine hands over up to ioBatchSize records per
	// channel send instead of one, so this hot path pays one channel
	// send/receive per batch rather than one per record.
	recChan, br := xam.NewXamReaderChanBatch(infile, xam.ChanCap, xam.IoBuf, 1, ioBatchSize)

	// ── establish temp dir -------------───────────────────────────────────────────────
	if tmpDir != "" {
		if err := os.MkdirAll(tmpDir, 0o755); err != nil {
			return err
		}
	} else {
		tmpDir = "."
	}

	// ── update PG line in header ───────────────────────────────────────────────
	hdr := br.Header()
	pg := xam.NewPG(hdr, "djinn sam sort "+infile)
	if err := hdr.AddProgram(pg); err != nil {
		return err
	}

	chunks := make(chan []*sam.Record, threads)
	results := make(chan string, threads)
	errCh := make(chan error, threads)

	var wg sync.WaitGroup
	for range threads {
		wg.Go(func() {
			for c := range chunks {
				path, err := sortAndSpill(c, hdr, tmpDir)
				if err != nil {
					select {
					case errCh <- err:
					default:
					}
					continue
				}
				results <- path
			}
		})
	}

	// collector goroutine: drain results as they arrive so workers never block
	var tmpFiles []string
	var collectWG sync.WaitGroup
	collectWG.Go(func() {
		for path := range results {
			tmpFiles = append(tmpFiles, path)
		}
	})

	// producer: read + batch on the calling goroutine
	var chunk []*sam.Record

	for batch := range recChan {
		chunk = append(chunk, batch...)
		if len(chunk) >= chunkSize {
			chunks <- chunk
			chunk = nil
		}
	}

	if len(chunk) > 0 {
		chunks <- chunk
	}
	close(chunks)

	wg.Wait()
	close(results)
	collectWG.Wait()

	// cleanup temp files on any exit path
	defer func() {
		for _, f := range tmpFiles {
			os.Remove(f)
		}
		os.Remove(tmpDir)
	}()

	select {
	case err := <-errCh:
		return err
	default:
	}

	return mergeChunks(tmpFiles, hdr, outPath, threads, asSam)
}

// bxTagName is the two-byte SAM tag ID looked up on every record; hoisted to
// a var so the lookup never allocates a fresh []byte literal per call.
var bxTagName = []byte("BX")

// bxKeyOf returns r's BX:Z value as a zero-copy view into the record's own
// AuxFields backing array -- no string conversion, no interface boxing.
// sam.Aux.Value() would do both (string(a[3:]) allocates and copies, and
// the interface{} return boxes it again); since all every caller here does
// with the result is bytes.Compare it, there's nothing to buy by going
// through either step.
//
// Call this once per record, not once per comparison: Record.Tag does a
// linear scan of AuxFields, so calling it from inside a sort/heap
// comparator -- as this package used to -- means every record gets
// rescanned and reboxed on every comparison it takes part in (O(n log n)
// times during the sort, O(n log k) more times during the merge), not
// once. Both sortAndSpill and mergeChunks below compute each record's key
// exactly once and carry it alongside the record from then on.
func bxKeyOf(r *sam.Record) []byte {
	if aux, ok := r.Tag(bxTagName); ok {
		return aux[3:] // raw 'Z'-type value bytes, same slice Value() would copy from
	}
	return missingBXKey
}

// keyedRecord pairs a record with its precomputed sort key.
type keyedRecord struct {
	key []byte
	rec *sam.Record
}

// sortAndSpill sorts one chunk by BX:Z and writes it to a fresh temp BAM. Each
// call owns its own slice, so no lock is needed across worker goroutines.
func sortAndSpill(chunk []*sam.Record, hdr *sam.Header, tmpDir string) (string, error) {
	keyed := make([]keyedRecord, len(chunk))
	for i, r := range chunk {
		keyed[i] = keyedRecord{key: bxKeyOf(r), rec: r}
	}
	// slices.SortFunc on a concrete []keyedRecord avoids sort.Slice's
	// reflection-based swap, and the comparator is now a plain
	// bytes.Compare against an already-computed key -- no further tag
	// lookups.
	slices.SortFunc(keyed, func(a, b keyedRecord) int { return bytes.Compare(a.key, b.key) })

	f, err := os.CreateTemp(tmpDir, ".bxsort-*.bam")
	if err != nil {
		return "", err
	}
	defer f.Close()

	bw, err := bam.NewWriterLevel(f, hdr, 0, 0)
	if err != nil {
		return "", err
	}
	for _, kr := range keyed {
		if err := bw.Write(kr.rec); err != nil {
			bw.Close()
			return "", err
		}
	}
	if err := bw.Close(); err != nil {
		return "", err
	}
	return f.Name(), nil
}

// mergeItem is one entry in the k-way merge heap: a record, its
// precomputed sort key, and which temp file it came from.
type mergeItem struct {
	rec *sam.Record
	key []byte
	src int
}

// mergeHeap is a hand-rolled binary min-heap over a concrete []mergeItem.
// container/heap's Push/Pop take an `any`, so every record crossing it
// would be boxed into an interface and type-asserted back out; with up to
// a few hundred thousand records flowing through this heap per merge,
// that's a lot of avoidable allocation+assertion for no benefit, since the
// element type here is fixed and known at compile time.
type mergeHeap []mergeItem

func (h *mergeHeap) push(it mergeItem) {
	*h = append(*h, it)
	i := len(*h) - 1
	for i > 0 {
		parent := (i - 1) / 2
		if bytes.Compare((*h)[i].key, (*h)[parent].key) >= 0 {
			break
		}
		(*h)[i], (*h)[parent] = (*h)[parent], (*h)[i]
		i = parent
	}
}

func (h *mergeHeap) pop() mergeItem {
	old := *h
	n := len(old) - 1
	top := old[0]
	old[0] = old[n]
	*h = old[:n]
	h.siftDown(0)
	return top
}

func (h *mergeHeap) siftDown(i int) {
	n := len(*h)
	for {
		left := 2*i + 1
		if left >= n {
			return
		}
		smallest := left
		if right := left + 1; right < n && bytes.Compare((*h)[right].key, (*h)[left].key) < 0 {
			smallest = right
		}
		if bytes.Compare((*h)[smallest].key, (*h)[i].key) >= 0 {
			return
		}
		(*h)[i], (*h)[smallest] = (*h)[smallest], (*h)[i]
		i = smallest
	}
}

// mergeChunks does a single-threaded k-way merge of already-sorted temp
// files. This stage is I/O-bound with one writer, so it stays serial;
// parallelizing it would need a second merge pass over partial merges,
// more complexity than the win is worth.
func mergeChunks(tmpFiles []string, hdr *sam.Header, outPath string, threads int, asSam bool) (err error) {
	readers := make([]*bam.Reader, len(tmpFiles))
	files := make([]*os.File, len(tmpFiles))
	for i, path := range tmpFiles {
		f, ferr := os.Open(path)
		if ferr != nil {
			return ferr
		}
		files[i] = f
		br, berr := bam.NewReader(f, 0)
		if berr != nil {
			return berr
		}
		readers[i] = br
	}
	defer func() {
		for _, f := range files {
			f.Close()
		}
	}()

	// --- Create writer channel --------------------------
	// outPath is passed through verbatim ("-" means stdout); NewXamWriterChanBatch
	// owns opening/creating the destination file itself. Batched: records popped
	// off the heap are accumulated into a slice and sent ioBatchSize at a time,
	// rather than one channel send per record.
	writeChan, writeDone := xam.NewXamWriterChanBatch(outPath, hdr, xam.ChanCap, xam.IoBuf, threads, asSam)

	h := &mergeHeap{}
	for i, br := range readers {
		rec, rerr := br.Read()
		if rerr == nil {
			h.push(mergeItem{rec: rec, key: bxKeyOf(rec), src: i})
		} else if !errors.Is(rerr, io.EOF) {
			return rerr
		}
	}

	batch := make([]*sam.Record, 0, ioBatchSize)
	for len(*h) > 0 {
		top := h.pop()
		batch = append(batch, top.rec)
		if len(batch) >= ioBatchSize {
			writeChan <- batch
			batch = make([]*sam.Record, 0, ioBatchSize)
		}
		rec, rerr := readers[top.src].Read()
		if rerr == nil {
			h.push(mergeItem{rec: rec, key: bxKeyOf(rec), src: top.src})
		} else if !errors.Is(rerr, io.EOF) {
			return rerr
		}
	}
	if len(batch) > 0 {
		writeChan <- batch
	}
	close(writeChan)
	<-writeDone

	return nil
}
