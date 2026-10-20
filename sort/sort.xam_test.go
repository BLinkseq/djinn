package sort

import (
	"bytes"
	"io"
	mathrand "math/rand"
	"os"
	"sort"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ── helpers ──────────────────────────────────────────────────────────────────

func newTestHeader(t *testing.T) *sam.Header {
	t.Helper()
	hdr, err := sam.NewHeader(nil, nil)
	if err != nil {
		t.Fatalf("sam.NewHeader: %v", err)
	}
	return hdr
}

func newTestRecord(t *testing.T, name, bx string) *sam.Record {
	t.Helper()
	rec, err := sam.NewRecord(name, nil, nil, -1, -1, 0, 0, nil, []byte("ACGT"), []byte{40, 40, 40, 40}, nil)
	if err != nil {
		t.Fatalf("sam.NewRecord: %v", err)
	}
	rec.Flags = sam.Unmapped
	if bx != "" {
		aux, err := sam.NewAux(sam.Tag{'B', 'X'}, bx)
		if err != nil {
			t.Fatalf("sam.NewAux: %v", err)
		}
		rec.AuxFields = append(rec.AuxFields, aux)
	}
	return rec
}

func writeBamFile(t *testing.T, path string, hdr *sam.Header, recs []*sam.Record) {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("os.Create: %v", err)
	}
	defer f.Close()
	bw, err := bam.NewWriterLevel(f, hdr, 0, 1)
	if err != nil {
		t.Fatalf("bam.NewWriterLevel: %v", err)
	}
	for _, r := range recs {
		if err := bw.Write(r); err != nil {
			t.Fatalf("bw.Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bw.Close: %v", err)
	}
}

func readBamBX(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("os.Open: %v", err)
	}
	defer f.Close()
	br, err := bam.NewReader(f, 0)
	if err != nil {
		t.Fatalf("bam.NewReader: %v", err)
	}
	var out []string
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("br.Read: %v", err)
		}
		out = append(out, string(bxKeyOf(rec)))
	}
	return out
}

// ── bxKeyOf ──────────────────────────────────────────────────────────────────

func TestBxKeyOf(t *testing.T) {
	cases := []struct {
		name string
		bx   string
		want string
	}{
		{"has BX", "AAACCC", "AAACCC"},
		{"missing BX sorts last", "", "\xff\xff\xff\xff"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			rec := newTestRecord(t, "read1", c.bx)
			got := string(bxKeyOf(rec))
			if got != c.want {
				t.Errorf("bxKeyOf() = %q, want %q", got, c.want)
			}
		})
	}
}

func TestBxKeyOfOrdering(t *testing.T) {
	// A record with no BX tag must sort after any record with a real BX
	// value (samtools convention: unset tags trail).
	withBX := newTestRecord(t, "r1", "AAAA")
	withoutBX := newTestRecord(t, "r2", "")
	if !(bytes.Compare(bxKeyOf(withBX), bxKeyOf(withoutBX)) < 0) {
		t.Errorf("expected record with BX to sort before record without BX")
	}
}

// ── sortAndSpill ─────────────────────────────────────────────────────────────

func TestSortAndSpillOrdersByBX(t *testing.T) {
	hdr := newTestHeader(t)
	dir := t.TempDir()

	chunk := []*sam.Record{
		newTestRecord(t, "r1", "CCCC"),
		newTestRecord(t, "r2", "AAAA"),
		newTestRecord(t, "r3", ""), // missing BX, must land last
		newTestRecord(t, "r4", "BBBB"),
	}

	path, err := sortAndSpill(chunk, hdr, dir)
	if err != nil {
		t.Fatalf("sortAndSpill: %v", err)
	}
	defer os.Remove(path)

	got := readBamBX(t, path)
	want := []string{"AAAA", "BBBB", "CCCC", "\xff\xff\xff\xff"}
	if len(got) != len(want) {
		t.Fatalf("got %d records, want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("record %d: got %q, want %q (full: %v)", i, got[i], want[i], got)
		}
	}
}

func TestSortAndSpillEmptyChunk(t *testing.T) {
	hdr := newTestHeader(t)
	dir := t.TempDir()

	path, err := sortAndSpill(nil, hdr, dir)
	if err != nil {
		t.Fatalf("sortAndSpill with empty chunk should not error: %v", err)
	}
	defer os.Remove(path)

	got := readBamBX(t, path)
	if len(got) != 0 {
		t.Errorf("expected 0 records for empty chunk, got %d", len(got))
	}
}

// ── mergeHeap ────────────────────────────────────────────────────────────────

func mkMergeItem(t *testing.T, name, bx string, src int) mergeItem {
	t.Helper()
	rec := newTestRecord(t, name, bx)
	return mergeItem{rec: rec, key: bxKeyOf(rec), src: src}
}

func TestMergeHeapOrdering(t *testing.T) {
	// Pushing out of order must pop back out in ascending BX order.
	h := &mergeHeap{}
	h.push(mkMergeItem(t, "a", "BBBB", 0))
	h.push(mkMergeItem(t, "b", "AAAA", 1))
	h.push(mkMergeItem(t, "c", "CCCC", 2))

	var got []string
	for len(*h) > 0 {
		got = append(got, string(h.pop().key))
	}
	want := []string{"AAAA", "BBBB", "CCCC"}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("pop order = %v, want %v", got, want)
			break
		}
	}
}

func TestMergeHeapPushPop(t *testing.T) {
	h := &mergeHeap{}
	h.push(mkMergeItem(t, "x", "AAAA", 0))
	if len(*h) != 1 {
		t.Fatalf("len = %d, want 1", len(*h))
	}
	popped := h.pop()
	if string(popped.key) != "AAAA" {
		t.Errorf("pop() returned wrong item")
	}
	if len(*h) != 0 {
		t.Errorf("len after pop = %d, want 0", len(*h))
	}
}

// ── mergeChunks ──────────────────────────────────────────────────────────────

func TestMergeChunksWritesToOutPath(t *testing.T) {
	// Regression test: mergeChunks used to hardcode "-" (stdout) as the
	// writer destination regardless of the outPath argument, silently
	// discarding a real output path. This verifies outPath is honored.
	hdr := newTestHeader(t)
	dir := t.TempDir()

	chunkA := []*sam.Record{newTestRecord(t, "r1", "AAAA"), newTestRecord(t, "r3", "CCCC")}
	chunkB := []*sam.Record{newTestRecord(t, "r2", "BBBB")}

	pathA, err := sortAndSpill(chunkA, hdr, dir)
	if err != nil {
		t.Fatalf("sortAndSpill chunkA: %v", err)
	}
	pathB, err := sortAndSpill(chunkB, hdr, dir)
	if err != nil {
		t.Fatalf("sortAndSpill chunkB: %v", err)
	}

	outPath := dir + "/merged.bam"
	if err := mergeChunks([]string{pathA, pathB}, hdr, outPath, 2, false); err != nil {
		t.Fatalf("mergeChunks: %v", err)
	}

	info, err := os.Stat(outPath)
	if err != nil {
		t.Fatalf("expected output file to exist at outPath, but os.Stat failed: %v", err)
	}
	if info.Size() == 0 {
		t.Fatalf("output file at outPath is empty; merged records were not written there")
	}

	got := readBamBX(t, outPath)
	want := []string{"AAAA", "BBBB", "CCCC"}
	if len(got) != len(want) {
		t.Fatalf("got %d records, want %d (%v)", len(got), len(want), got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("record %d: got %q, want %q", i, got[i], want[i])
		}
	}
}

// ── SortByBX (end-to-end) ────────────────────────────────────────────────────

func TestSortByBXEndToEnd(t *testing.T) {
	dir := t.TempDir()
	hdr := newTestHeader(t)

	recs := []*sam.Record{
		newTestRecord(t, "r1", "TTTT"),
		newTestRecord(t, "r2", "AAAA"),
		newTestRecord(t, "r3", "GGGG"),
		newTestRecord(t, "r4", "CCCC"),
	}
	inPath := dir + "/in.bam"
	writeBamFile(t, inPath, hdr, recs)

	outPath := dir + "/out.bam"
	if err := SortByBX(inPath, outPath, dir, 2, false); err != nil {
		t.Fatalf("SortByBX: %v", err)
	}

	got := readBamBX(t, outPath)
	want := []string{"AAAA", "CCCC", "GGGG", "TTTT"}
	if len(got) != len(want) {
		t.Fatalf("got %d records, want %d (%v)", len(got), len(want), got)
	}
	if !sort.StringsAreSorted(got) {
		t.Errorf("output records are not sorted by BX: %v", got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("record %d: got %q, want %q", i, got[i], want[i])
		}
	}
}

func TestSortByBXEmptyInput(t *testing.T) {
	dir := t.TempDir()
	hdr := newTestHeader(t)

	inPath := dir + "/empty.bam"
	writeBamFile(t, inPath, hdr, nil)

	outPath := dir + "/out.bam"
	if err := SortByBX(inPath, outPath, dir, 2, false); err != nil {
		t.Fatalf("SortByBX on empty input should not error: %v", err)
	}

	got := readBamBX(t, outPath)
	if len(got) != 0 {
		t.Errorf("expected 0 records, got %d", len(got))
	}
}

// ── fuzz-ish correctness check for the hand-rolled heap ────────────────────

func TestMergeHeapRandomizedOrdering(t *testing.T) {
	r := mathrand.New(mathrand.NewSource(1))
	const n = 2000
	h := &mergeHeap{}
	var keys []string
	for i := 0; i < n; i++ {
		bx := make([]byte, 8)
		for j := range bx {
			bx[j] = byte('A' + r.Intn(26))
		}
		keys = append(keys, string(bx))
		h.push(mkMergeItem(t, "r", string(bx), i))
	}
	sort.Strings(keys)

	for i := 0; i < n; i++ {
		got := string(h.pop().key)
		if got != keys[i] {
			t.Fatalf("pop %d = %q, want %q", i, got, keys[i])
		}
	}
	if len(*h) != 0 {
		t.Errorf("heap not empty after popping all items: %d left", len(*h))
	}
}
