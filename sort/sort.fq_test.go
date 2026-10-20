package sort

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"github.com/blinkseq/djinn/fastq"
	"io"
	mathrand "math/rand"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// ── helpers ──────────────────────────────────────────────────────────────────

type hapRead struct {
	id string
	bc string // "" -> no BX tag on either mate
}

// writeHapFastqPair writes matching haplotagging-format R1/R2 FASTQ files:
// same read ID on both mates, same BX:Z tag (when bc != "") on both.
func writeHapFastqPair(t *testing.T, r1Path, r2Path string, reads []hapRead) {
	t.Helper()
	var r1, r2 strings.Builder
	for i, r := range reads {
		id := r.id
		if id == "" {
			id = fmt.Sprintf("read%d", i)
		}
		tag := ""
		if r.bc != "" {
			tag = " BX:Z:" + r.bc
		}
		fmt.Fprintf(&r1, "@%s%s\nACGTACGTAC\n+\nIIIIIIIIII\n", id, tag)
		fmt.Fprintf(&r2, "@%s%s\nTGCATGCATG\n+\nIIIIIIIIII\n", id, tag)
	}
	if err := os.WriteFile(r1Path, []byte(r1.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(r2Path, []byte(r2.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

// writeHapFastqSingle writes a single haplotagging-format FASTQ file (no
// mate), for single-end sort tests.
func writeHapFastqSingle(t *testing.T, path string, reads []hapRead) {
	t.Helper()
	var sb strings.Builder
	for i, r := range reads {
		id := r.id
		if id == "" {
			id = fmt.Sprintf("read%d", i)
		}
		tag := ""
		if r.bc != "" {
			tag = " BX:Z:" + r.bc
		}
		fmt.Fprintf(&sb, "@%s%s\nACGTACGTAC\n+\nIIIIIIIIII\n", id, tag)
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

// readGzFastqIDsAndBX reads a gzip FASTQ file and returns, in record order,
// each record's ID and its BX:Z value (empty string if none).
func readGzFastqIDsAndBX(t *testing.T, path string) (ids, bxs []string) {
	t.Helper()
	fh, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	gz, err := gzip.NewReader(fh)
	if err != nil {
		t.Fatal(err)
	}
	defer gz.Close()

	s := bufio.NewScanner(gz)
	s.Buffer(make([]byte, 1<<20), 1<<20)
	for s.Scan() {
		line := s.Text()
		if !strings.HasPrefix(line, "@") {
			continue
		}
		header := strings.TrimPrefix(line, "@")
		fields := strings.Fields(header)
		// Format-preserving (non-standardize) output appends "/1" or "/2"
		// to the ID (ToHaplotagging's mate-digit suffix); strip it so
		// assertions stay agnostic to which converter wrote the record.
		id := strings.TrimSuffix(strings.TrimSuffix(fields[0], "/1"), "/2")
		ids = append(ids, id)
		bx := ""
		for _, f := range fields[1:] {
			if after, ok := strings.CutPrefix(f, "BX:Z:"); ok {
				bx = after
			}
		}
		bxs = append(bxs, bx)
		// skip seq, '+', qual
		for range 3 {
			if !s.Scan() {
				t.Fatalf("truncated FASTQ record after %q", line)
			}
		}
		if s.Err() != nil {

		}
	}
	return ids, bxs
}

// ── end-to-end ───────────────────────────────────────────────────────────────

func TestSortFqByBX_SortsByBarcodeAndKeepsPairsMatched(t *testing.T) {
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")

	reads := []hapRead{
		{id: "c", bc: "CCCC"},
		{id: "a", bc: "AAAA"},
		{id: "none", bc: ""}, // no barcode, must sort last
		{id: "b", bc: "BBBB"},
	}
	writeHapFastqPair(t, r1, r2, reads)

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, ""); err != nil {
		t.Fatalf("SortFqByBX: %v", err)
	}

	gotIDs1, gotBX1 := readGzFastqIDsAndBX(t, prefix+".R1.fq.gz")
	gotIDs2, gotBX2 := readGzFastqIDsAndBX(t, prefix+".R2.fq.gz")

	wantIDs := []string{"a", "b", "c", "none"}
	wantBX := []string{"AAAA", "BBBB", "CCCC", ""}

	if len(gotIDs1) != len(wantIDs) {
		t.Fatalf("R1: got %d records, want %d (%v)", len(gotIDs1), len(wantIDs), gotIDs1)
	}
	for i := range wantIDs {
		if gotIDs1[i] != wantIDs[i] || gotBX1[i] != wantBX[i] {
			t.Errorf("R1 record %d: got (id=%q,bx=%q), want (id=%q,bx=%q)", i, gotIDs1[i], gotBX1[i], wantIDs[i], wantBX[i])
		}
	}
	// R2 must be in the exact same order as R1 -- pairs travel together
	// through the sort, not independently.
	if !equalSlices(gotIDs1, gotIDs2) {
		t.Errorf("R1/R2 ID order mismatch: R1=%v R2=%v", gotIDs1, gotIDs2)
	}
	if !equalSlices(gotBX1, gotBX2) {
		t.Errorf("R1/R2 BX order mismatch: R1=%v R2=%v", gotBX1, gotBX2)
	}
}

func TestSortFqByBX_SingleEndSortsByBarcode(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")

	reads := []hapRead{
		{id: "c", bc: "CCCC"},
		{id: "a", bc: "AAAA"},
		{id: "none", bc: ""}, // no barcode, must sort last
		{id: "b", bc: "BBBB"},
	}
	writeHapFastqSingle(t, r1, reads)

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1}, prefix, dir, 2, false, ""); err != nil {
		t.Fatalf("SortFqByBX: %v", err)
	}

	gotIDs, gotBX := readGzFastqIDsAndBX(t, prefix+".R1.fq.gz")
	wantIDs := []string{"a", "b", "c", "none"}
	wantBX := []string{"AAAA", "BBBB", "CCCC", ""}

	if len(gotIDs) != len(wantIDs) {
		t.Fatalf("got %d records, want %d (%v)", len(gotIDs), len(wantIDs), gotIDs)
	}
	for i := range wantIDs {
		if gotIDs[i] != wantIDs[i] || gotBX[i] != wantBX[i] {
			t.Errorf("record %d: got (id=%q,bx=%q), want (id=%q,bx=%q)", i, gotIDs[i], gotBX[i], wantIDs[i], wantBX[i])
		}
	}
	if _, err := os.Stat(prefix + ".R2.fq.gz"); err == nil {
		t.Errorf("expected no R2 output file for single-end input")
	}
}

func TestSortFqByBX_SingleEndDoesNotForceCasavaMarker(t *testing.T) {
	// Unlike paired output, single-end output has no "other side" to
	// distinguish from, so the mate digit must never be synthesized --
	// only written if the source record actually carried one.
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	writeHapFastqSingle(t, r1, []hapRead{{id: "a", bc: "AAAA"}})

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1}, prefix, dir, 2, false, ""); err != nil {
		t.Fatalf("SortFqByBX: %v", err)
	}

	raw := readGzFirstHeaderLine(t, prefix+".R1.fq.gz")
	if raw != "@a\tBX:Z:AAAA" {
		t.Errorf("expected no synthesized CASAVA suffix on single-end output, got %q", raw)
	}
}

func TestSortFqByBX_DefaultPreservesInputFormat(t *testing.T) {
	// standardize=false (the default): haplotagging input should come
	// back out via ToHaplotagging -- a mate digit is always appended to
	// the ID ("/1"/"/2"), even though none was present on the source
	// records, since that's how haplotagging format represents mate
	// pairing and this package's bare BX tags carry no CASAVA marker.
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")
	writeHapFastqPair(t, r1, r2, []hapRead{{id: "a", bc: "AAAA"}})

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, ""); err != nil {
		t.Fatalf("SortFqByBX: %v", err)
	}

	raw := readGzFirstHeaderLine(t, prefix+".R1.fq.gz")
	if raw != "@a/1\tBX:Z:AAAA" {
		t.Errorf("expected ToHaplotagging-style header, got %q", raw)
	}
	raw2 := readGzFirstHeaderLine(t, prefix+".R2.fq.gz")
	if raw2 != "@a/2\tBX:Z:AAAA" {
		t.Errorf("expected ToHaplotagging-style header, got %q", raw2)
	}
}

func TestSortFqByBX_StandardizeFlag(t *testing.T) {
	// standardize=true: output uses ToStandard's generic BX:Z:/VX:i: form,
	// and an absent CASAVA marker stays absent (no synthesized "/1"/"/2").
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")
	writeHapFastqPair(t, r1, r2, []hapRead{{id: "a", bc: "AAAA"}})

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, true, ""); err != nil {
		t.Fatalf("SortFqByBX: %v", err)
	}

	raw := readGzFirstHeaderLine(t, prefix+".R1.fq.gz")
	if raw != "@a\tBX:Z:AAAA\tVX:i:1" {
		t.Errorf("expected ToStandard-style header, got %q", raw)
	}
}

func TestSortFqByBX_PreservesFullCasavaTextForStlfr(t *testing.T) {
	// ToStlfr (and ToTellseq/ToTenX) write a mate's entire CASAVA marker
	// back out verbatim when it's a new-style Illumina field (longer
	// than the bare '1'/'2' digit ToStandard/ToHaplotagging settle for).
	// A record's real filter flag and index sequence must survive the
	// sort unchanged -- not collapse to a generic reconstructed stand-in.
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")
	if err := os.WriteFile(r1, []byte("@read1#123_456_789 1:N:0:ATCGATCGATCG\nACGTACGTAC\n+\nIIIIIIIIII\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(r2, []byte("@read1#123_456_789 2:N:0:ATCGATCGATCG\nTGCATGCATG\n+\nIIIIIIIIII\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, ""); err != nil {
		t.Fatalf("SortFqByBX: %v", err)
	}

	raw1 := readGzFirstHeaderLine(t, prefix+".R1.fq.gz")
	if want := "@read1#123_456_789\t1:N:0:ATCGATCGATCG"; raw1 != want {
		t.Errorf("R1 header = %q, want %q (original CASAVA marker preserved verbatim)", raw1, want)
	}
	raw2 := readGzFirstHeaderLine(t, prefix+".R2.fq.gz")
	if want := "@read1#123_456_789\t2:N:0:ATCGATCGATCG"; raw2 != want {
		t.Errorf("R2 header = %q, want %q (original CASAVA marker preserved verbatim)", raw2, want)
	}
}

func TestSortFqByBX_FallsBackTo10XWhenBcMapProvided(t *testing.T) {
	// No BX:Z:/stlfr/tellseq marker is present on these records, so
	// CoreFqParser can't detect a format; with a 10X whitelist provided,
	// SortFqByBX should fall back to the 10X path (barcode inline as a
	// 16nt R1 sequence prefix), same dispatch convert.ConvertFq uses.
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")

	const bc = "AAAACCCCGGGGTTTT" // 16nt
	os.WriteFile(r1, []byte("@read1\n"+bc+"ACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"), 0o644)
	os.WriteFile(r2, []byte("@read1\nTGCA\n+\nJJJJ\n"), 0o644)

	bcMap := filepath.Join(dir, "whitelist.txt")
	if err := os.WriteFile(bcMap, []byte(bc+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, bcMap); err != nil {
		t.Fatalf("SortFqByBX: %v", err)
	}

	raw1 := readGzFirstHeaderLine(t, prefix+".R1.fq.gz")
	if raw1 != "@read1:"+bc+"\t1:N:18:ATCG" {
		t.Errorf("expected ToTenX-style R1 header, got %q", raw1)
	}
	raw2 := readGzFirstHeaderLine(t, prefix+".R2.fq.gz")
	if raw2 != "@read1:"+bc+"\t2:N:18:ATCG" {
		t.Errorf("expected ToTenX-style R2 header, got %q", raw2)
	}
}

func TestSortFqByBX_NoFormatDetectedAndNoBcMapIsError(t *testing.T) {
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")
	os.WriteFile(r1, []byte("@read1\nACGTACGTAC\n+\nIIIIIIIIII\n"), 0o644)
	os.WriteFile(r2, []byte("@read1\nTGCATGCATG\n+\nIIIIIIIIII\n"), 0o644)

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, ""); err == nil {
		t.Fatal("expected an error when no format is detected and no bc-map is given")
	}
}

func readGzFirstHeaderLine(t *testing.T, path string) string {
	t.Helper()
	fh, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	gz, err := gzip.NewReader(fh)
	if err != nil {
		t.Fatal(err)
	}
	defer gz.Close()
	s := bufio.NewScanner(gz)
	if !s.Scan() {
		t.Fatal("empty output file")
	}
	return s.Text()
}

func TestSortFqByBX_MismatchedMateIDsIsError(t *testing.T) {
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")

	writeHapFastqPair(t, r1, r2, []hapRead{{id: "a", bc: "AAAA"}})
	// Overwrite r2 with a different ID for the same position.
	if err := os.WriteFile(r2, []byte("@different BX:Z:AAAA\nTGCATGCATG\n+\nIIIIIIIIII\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, ""); err == nil {
		t.Fatal("expected an error for mismatched mate IDs, got nil")
	}
}

func TestSortFqByBX_MismatchedMateCountIsError(t *testing.T) {
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")

	writeHapFastqPair(t, r1, r2, []hapRead{{id: "a", bc: "AAAA"}, {id: "b", bc: "BBBB"}})
	// Truncate r2 to only one record.
	if err := os.WriteFile(r2, []byte("@a BX:Z:AAAA\nTGCATGCATG\n+\nIIIIIIIIII\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	prefix := filepath.Join(dir, "out")
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, ""); err == nil {
		t.Fatal("expected an error for mismatched mate counts, got nil")
	}
}

func TestSortFqByBX_EmptyInput(t *testing.T) {
	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "r1.fq"), filepath.Join(dir, "r2.fq")
	writeHapFastqPair(t, r1, r2, nil)

	prefix := filepath.Join(dir, "out")
	// With zero records, no format can be sampled/detected -- same
	// behavior as every other fastq.CoreFqParser caller (ConvertFq,
	// FilterInvalidFQ, ...) on a fully empty input.
	if err := SortFqByBX([]string{r1, r2}, prefix, dir, 2, false, ""); err == nil {
		t.Fatal("expected an error when no format can be detected from empty input, got nil")
	}
}

func TestMergeFqChunks_ManySpilledChunks(t *testing.T) {
	// Exercise the k-way merge path directly (several small spilled
	// chunks) rather than forcing SortFqByBX through a real chunkSize
	// (250,000-pair) boundary end-to-end, which is correct but far too
	// slow for a unit test.
	dir := t.TempDir()
	r := mathrand.New(mathrand.NewSource(1))

	const numChunks = 5
	const perChunk = 40
	var wantBX []string
	var tmpFiles []string
	for c := range numChunks {
		chunk := make([]fqPair, perChunk)
		for i := range chunk {
			bc := make([]byte, 6)
			for j := range bc {
				bc[j] = byte('A' + r.Intn(26))
			}
			chunk[i] = fqPair{
				id:      fmt.Appendf(nil, "c%d_r%d", c, i),
				barcode: append([]byte(nil), bc...),
				valid:   true,
				mates: []mate{
					{seq: []byte("ACGT"), qual: []byte("IIII")},
					{seq: []byte("TGCA"), qual: []byte("JJJJ")},
				},
			}
			wantBX = append(wantBX, string(bc))
		}
		path, err := fqsortAndSpill(chunk, dir)
		if err != nil {
			t.Fatalf("sortAndSpill: %v", err)
		}
		tmpFiles = append(tmpFiles, path)
	}
	sort.Strings(wantBX)

	prefix := filepath.Join(dir, "out")
	if err := mergeFqChunks(tmpFiles, prefix, (*fastq.CoreFq).ToStandard, true, 2); err != nil {
		t.Fatalf("mergeFqChunks: %v", err)
	}

	_, gotBX := readGzFastqIDsAndBX(t, prefix+".R1.fq.gz")
	_, gotBX2 := readGzFastqIDsAndBX(t, prefix+".R2.fq.gz")
	if len(gotBX) != len(wantBX) {
		t.Fatalf("got %d records, want %d", len(gotBX), len(wantBX))
	}
	for i := range wantBX {
		if gotBX[i] != wantBX[i] {
			t.Fatalf("R1 record %d: got BX %q, want %q", i, gotBX[i], wantBX[i])
		}
		if gotBX2[i] != wantBX[i] {
			t.Fatalf("R2 record %d: got BX %q, want %q", i, gotBX2[i], wantBX[i])
		}
	}
}

func equalSlices(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// ── binary spill format round-trip ──────────────────────────────────────────

func TestWriteReadPairRoundTrip(t *testing.T) {
	p := fqPair{
		id:       []byte("read1"),
		barcode:  []byte("AAAA"),
		comments: []byte("XX:i:5"),
		valid:    true,
		mates: []mate{
			{seq: []byte("ACGT"), qual: []byte("IIII"), casava: []byte("1")},
			{seq: []byte("TGCA"), qual: []byte("JJJJ"), casava: []byte("2")},
		},
	}
	var buf bytes.Buffer
	w := bufio.NewWriter(&buf)
	if err := writePair(w, &p); err != nil {
		t.Fatalf("writePair: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}

	got, err := readPair(bufio.NewReader(&buf), 2)
	if err != nil {
		t.Fatalf("readPair: %v", err)
	}
	if string(got.id) != "read1" || string(got.barcode) != "AAAA" || string(got.comments) != "XX:i:5" || !got.valid {
		t.Errorf("round trip mismatch: %+v", got)
	}
	if string(got.mates[0].seq) != "ACGT" || string(got.mates[0].qual) != "IIII" || string(got.mates[0].casava) != "1" {
		t.Errorf("mate1 round trip mismatch: %+v", got)
	}
	if string(got.mates[1].seq) != "TGCA" || string(got.mates[1].qual) != "JJJJ" || string(got.mates[1].casava) != "2" {
		t.Errorf("mate2 round trip mismatch: %+v", got)
	}
}

// TestWriteReadPairRoundTrip_FullCasavaTextPreserved guards the bug this
// struct's earlier hasCasava-bool design had: ToStlfr/ToTellseq/ToTenX
// write a mate's *entire* CASAVA marker back out when it's longer than
// one byte, not just its first byte like ToStandard/ToHaplotagging do.
// A new-style Illumina marker's real filter flag and index sequence
// must survive the spill/merge round trip verbatim, not collapse to a
// bare "1"/"2" and get reconstructed as a generic placeholder on write.
func TestWriteReadPairRoundTrip_FullCasavaTextPreserved(t *testing.T) {
	p := fqPair{
		id:      []byte("read1"),
		barcode: []byte("AAAA"),
		mates: []mate{
			{seq: []byte("ACGT"), qual: []byte("IIII"), casava: []byte("1:N:0:ATCGATCGATCG")},
			{seq: []byte("TGCA"), qual: []byte("JJJJ"), casava: []byte("2:N:0:ATCGATCGATCG")},
		},
	}
	var buf bytes.Buffer
	w := bufio.NewWriter(&buf)
	if err := writePair(w, &p); err != nil {
		t.Fatalf("writePair: %v", err)
	}
	w.Flush()

	got, err := readPair(bufio.NewReader(&buf), 2)
	if err != nil {
		t.Fatalf("readPair: %v", err)
	}
	if string(got.mates[0].casava) != "1:N:0:ATCGATCGATCG" {
		t.Errorf("mate1 CASAVA = %q, want the original marker verbatim", got.mates[0].casava)
	}
	if string(got.mates[1].casava) != "2:N:0:ATCGATCGATCG" {
		t.Errorf("mate2 CASAVA = %q, want the original marker verbatim", got.mates[1].casava)
	}
}

func TestWriteReadPairRoundTrip_NoCasavaOrComments(t *testing.T) {
	p := fqPair{
		id:      []byte("read1"),
		barcode: []byte("AAAA"),
		mates: []mate{
			{seq: []byte("ACGT"), qual: []byte("IIII")},
			{seq: []byte("TGCA"), qual: []byte("JJJJ")},
		},
	}
	var buf bytes.Buffer
	w := bufio.NewWriter(&buf)
	if err := writePair(w, &p); err != nil {
		t.Fatalf("writePair: %v", err)
	}
	w.Flush()

	got, err := readPair(bufio.NewReader(&buf), 2)
	if err != nil {
		t.Fatalf("readPair: %v", err)
	}
	if len(got.mates[0].casava) != 0 || len(got.mates[1].casava) != 0 {
		t.Errorf("expected no CASAVA marker on either mate, got %+v", got)
	}
	if len(got.comments) != 0 {
		t.Errorf("expected empty comments, got %q", got.comments)
	}
}

func TestWriteReadPairRoundTrip_SingleEnd(t *testing.T) {
	p := fqPair{
		id:       []byte("read1"),
		barcode:  []byte("AAAA"),
		comments: []byte("XX:i:5"),
		valid:    true,
		mates: []mate{
			{seq: []byte("ACGT"), qual: []byte("IIII"), casava: []byte("1")},
		},
	}
	var buf bytes.Buffer
	w := bufio.NewWriter(&buf)
	if err := writePair(w, &p); err != nil {
		t.Fatalf("writePair: %v", err)
	}
	w.Flush()

	got, err := readPair(bufio.NewReader(&buf), 1)
	if err != nil {
		t.Fatalf("readPair: %v", err)
	}
	if len(got.mates) != 1 {
		t.Fatalf("got %d mates, want 1", len(got.mates))
	}
	if string(got.mates[0].seq) != "ACGT" || string(got.mates[0].qual) != "IIII" || string(got.mates[0].casava) != "1" {
		t.Errorf("mate round trip mismatch: %+v", got)
	}
}

func TestReadPairEOFAtBoundary(t *testing.T) {
	var buf bytes.Buffer
	w := bufio.NewWriter(&buf)
	p := fqPair{id: []byte("a"), barcode: []byte("AAAA")}
	if err := writePair(w, &p); err != nil {
		t.Fatal(err)
	}
	w.Flush()

	r := bufio.NewReader(&buf)
	if _, err := readPair(r, 0); err != nil {
		t.Fatalf("first readPair: %v", err)
	}
	if _, err := readPair(r, 0); err != io.EOF {
		t.Fatalf("second readPair: got %v, want io.EOF", err)
	}
}

// ── fqMergeHeap ──────────────────────────────────────────────────────────────

func TestFqMergeHeapOrdering(t *testing.T) {
	h := &fqMergeHeap{}
	h.push(fqMergeItem{pair: &fqPair{barcode: []byte("BBBB")}, src: 0})
	h.push(fqMergeItem{pair: &fqPair{barcode: []byte("AAAA")}, src: 1})
	h.push(fqMergeItem{pair: &fqPair{barcode: []byte("CCCC")}, src: 2})

	var got []string
	for len(*h) > 0 {
		got = append(got, string(h.pop().pair.barcode))
	}
	want := []string{"AAAA", "BBBB", "CCCC"}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("pop order = %v, want %v", got, want)
			break
		}
	}
}

func TestFqMergeHeapMissingBarcodeSortsLast(t *testing.T) {
	h := &fqMergeHeap{}
	h.push(fqMergeItem{pair: &fqPair{barcode: nil}, src: 0})
	h.push(fqMergeItem{pair: &fqPair{barcode: []byte("AAAA")}, src: 1})

	first := h.pop()
	if string(first.pair.barcode) != "AAAA" {
		t.Errorf("expected barcoded pair to pop first, got %+v", first.pair)
	}
}
