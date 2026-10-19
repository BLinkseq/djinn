package count

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/biogo/hts/sam"
)

// hapRead describes one synthetic haplotagging-format FASTQ record.
type hapRead struct {
	id  string
	bc  string // "" -> no BX tag
	vx  int    // -1 -> no VX tag, else 0/1
	seq string
}

// writeHapFastq writes a plain-text (uncompressed) FASTQ file containing the
// given records, in Haplotagging BX:Z:/VX:i: format.
func writeHapFastq(t *testing.T, path string, reads []hapRead) {
	t.Helper()
	var sb strings.Builder
	for i, r := range reads {
		seq := r.seq
		if seq == "" {
			seq = "ACGTACGTAC"
		}
		qual := strings.Repeat("I", len(seq))
		var desc strings.Builder
		if r.bc != "" {
			desc.WriteString(" BX:Z:")
			desc.WriteString(r.bc)
		}
		if r.vx >= 0 {
			fmt.Fprintf(&desc, " VX:i:%d", r.vx)
		}
		id := r.id
		if id == "" {
			id = fmt.Sprintf("read%d", i)
		}
		sb.WriteString("@" + id + desc.String() + "\n")
		sb.WriteString(seq + "\n+\n" + qual + "\n")
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

// padHapReads pads reads with harmless, uniquely-barcoded filler records so
// the file has at least `total` records. This is required because
// fastq.CoreFqParser insists on reading (and succeeding on) the first
// 100 records of the file to sniff the linked-read technology.
func padHapReads(reads []hapRead, total int) []hapRead {
	for i := 0; len(reads) < total; i++ {
		reads = append(reads, hapRead{bc: fmt.Sprintf("FILLER%04d", i), vx: 1})
	}
	return reads
}

// captureStdout redirects os.Stdout for the duration of fn and returns
// whatever was written to it, along with fn's error.
func captureStdout(t *testing.T, fn func() error) (string, error) {
	t.Helper()
	orig := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	os.Stdout = w
	outCh := make(chan string, 1)
	go func() {
		var buf bytes.Buffer
		io.Copy(&buf, r)
		outCh <- buf.String()
	}()

	fnErr := fn()

	w.Close()
	os.Stdout = orig
	out := <-outCh
	return out, fnErr
}

// parseCounts parses the tab-separated "barcode\tcount" output produced by
// CountFQ/CountXam into a map.
func parseCounts(t *testing.T, out string) map[string]int {
	t.Helper()
	counts := make(map[string]int)
	for _, line := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) != 2 {
			t.Fatalf("unexpected output line %q", line)
		}
		var v int
		if _, err := fmt.Sscanf(fields[1], "%d", &v); err != nil {
			t.Fatalf("bad count in line %q: %v", line, err)
		}
		counts[fields[0]] = v
	}
	return counts
}

// samRec is a minimal description of a SAM alignment record for test fixtures.
type samRec struct {
	name string
	bx   string // "" -> no BX tag
	vx   int    // -1 -> no VX tag, else 0/1
}

// writeSamFile writes a minimal, valid, unmapped-reads SAM text file
// containing the given records.
func writeSamFile(t *testing.T, path string, recs []samRec) {
	t.Helper()
	var sb strings.Builder
	sb.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	sb.WriteString("@SQ\tSN:chr1\tLN:1000\n")
	for _, r := range recs {
		var tags strings.Builder
		if r.bx != "" {
			tags.WriteString("\tBX:Z:")
			tags.WriteString(r.bx)
		}
		if r.vx >= 0 {
			fmt.Fprintf(&tags, "\tVX:i:%d", r.vx)
		}
		sb.WriteString(fmt.Sprintf("%s\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII%s\n", r.name, tags.String()))
	}
	if err := os.WriteFile(path, []byte(sb.String()), 0o644); err != nil {
		t.Fatal(err)
	}
}

// readAllSamRecords opens and fully reads a SAM/BAM file, returning all
// records. Used to verify writer output.
func readAllSamRecords(t *testing.T, path string) []*sam.Record {
	t.Helper()
	fh, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	r, err := sam.NewReader(fh)
	if err != nil {
		t.Fatal(err)
	}
	var out []*sam.Record
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, rec)
	}
	return out
}

// tempPath returns a path inside t.TempDir() with the given filename.
func tempPath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join(t.TempDir(), name)
}
