package invalid

import (
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/biogo/hts/sam"
)

type hapRead struct {
	id  string
	bc  string
	vx  int
	seq string
}

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

func padHapReads(reads []hapRead, total int) []hapRead {
	for i := 0; len(reads) < total; i++ {
		reads = append(reads, hapRead{bc: fmt.Sprintf("FILLER%04d", i), vx: 1})
	}
	return reads
}

// readGzFastqIDs reads a gzip-compressed FASTQ file and returns the record
// IDs it contains (without the leading '@' or description).
func readGzFastqIDs(t *testing.T, path string) []string {
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
	data, err := io.ReadAll(gz)
	if err != nil {
		t.Fatal(err)
	}
	var ids []string
	lines := strings.Split(string(data), "\n")
	for i := 0; i < len(lines); i += 4 {
		if lines[i] == "" {
			continue
		}
		header := strings.TrimPrefix(lines[i], "@")
		ids = append(ids, strings.Fields(header)[0])
	}
	return ids
}

type samRec struct {
	name string
	bx   string
	vx   int
}

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

func recNames(recs []*sam.Record) []string {
	names := make([]string, len(recs))
	for i, r := range recs {
		names[i] = r.Name
	}
	return names
}

func tempPath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join(t.TempDir(), name)
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
		data, _ := io.ReadAll(r)
		outCh <- string(data)
	}()

	fnErr := fn()

	w.Close()
	os.Stdout = orig
	out := <-outCh
	return out, fnErr
}

func parseSamText(t *testing.T, text string) []*sam.Record {
	t.Helper()
	r, err := sam.NewReader(strings.NewReader(text))
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
