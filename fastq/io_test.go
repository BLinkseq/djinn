package fastq

import (
	"bufio"
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/klauspost/pgzip"
)

func readGzipFile(t *testing.T, path string) []byte {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("opening %s: %v", path, err)
	}
	defer f.Close()
	gr, err := pgzip.NewReader(f)
	if err != nil {
		t.Fatalf("creating gzip reader: %v", err)
	}
	defer gr.Close()
	data, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("reading gzip content: %v", err)
	}
	return data
}

func TestFastqWriter_WriteRecordNoAux(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.fq.gz")
	fw, err := NewFastqWriter(path, "/1", 1, 1)
	if err != nil {
		t.Fatalf("NewFastqWriter: %v", err)
	}

	seq := []byte("ACGTACGTAC")
	qual := []uint8{40, 40, 40, 40, 40, 40, 40, 40, 40, 40} // Phred 40 -> 'I'
	if err := fw.WriteRecord("read1", nil, seq, qual); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := fw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	got := readGzipFile(t, path)
	want := "@read1/1\nACGTACGTAC\n+\nIIIIIIIIII\n"
	if string(got) != want {
		t.Errorf("record = %q, want %q", string(got), want)
	}
}

func TestFastqWriter_WriteRecordWithAux(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.fq.gz")
	fw, err := NewFastqWriter(path, "/2", 1, 1)
	if err != nil {
		t.Fatalf("NewFastqWriter: %v", err)
	}

	bxAux, err := sam.NewAux(sam.NewTag("BX"), "SAMPLEBC")
	if err != nil {
		t.Fatalf("sam.NewAux: %v", err)
	}
	aux := sam.AuxFields{bxAux}

	seq := []byte("ACGT")
	qual := []uint8{0, 1, 2, 3}
	if err := fw.WriteRecord("read2", aux, seq, qual); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := fw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	got := readGzipFile(t, path)

	// The header line must have a tab separating the read name/dir from the
	// first (and here, only) aux field -- this is the fix for the bug where
	// no separator was written before the first aux field.
	nl := bytes.IndexByte(got, '\n')
	if nl < 0 {
		t.Fatalf("no newline found in output: %q", got)
	}
	header := string(got[:nl])
	wantHeader := "@read2/2\t" + bxAux.String()
	if header != wantHeader {
		t.Errorf("header = %q, want %q", header, wantHeader)
	}

	wantQual := []byte{0 + 33, 1 + 33, 2 + 33, 3 + 33}
	wantBody := header + "\nACGT\n+\n" + string(wantQual) + "\n"
	if string(got) != wantBody {
		t.Errorf("record = %q, want %q", string(got), wantBody)
	}
}

func TestFastqWriter_WriteRecordMultipleAuxFields(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.fq.gz")
	fw, err := NewFastqWriter(path, "/1", 1, 1)
	if err != nil {
		t.Fatalf("NewFastqWriter: %v", err)
	}

	bxAux, err := sam.NewAux(sam.NewTag("BX"), "SAMPLEBC")
	if err != nil {
		t.Fatalf("sam.NewAux BX: %v", err)
	}
	vxAux, err := sam.NewAux(sam.NewTag("VX"), 1)
	if err != nil {
		t.Fatalf("sam.NewAux VX: %v", err)
	}
	aux := sam.AuxFields{vxAux, bxAux}

	if err := fw.WriteRecord("read3", aux, []byte("AC"), []uint8{5, 5}); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := fw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	got := readGzipFile(t, path)
	nl := bytes.IndexByte(got, '\n')
	header := string(got[:nl])
	wantHeader := "@read3/1\t" + vxAux.String() + "\t" + bxAux.String()
	if header != wantHeader {
		t.Errorf("header = %q, want %q", header, wantHeader)
	}
}

func TestFastqWriter_MultipleRecords(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.fq.gz")
	fw, err := NewFastqWriter(path, "/1", 1, 1)
	if err != nil {
		t.Fatalf("NewFastqWriter: %v", err)
	}
	for i := 0; i < 3; i++ {
		if err := fw.WriteRecord("read", nil, []byte("AC"), []uint8{1, 1}); err != nil {
			t.Fatalf("WriteRecord %d: %v", i, err)
		}
	}
	if err := fw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	got := readGzipFile(t, path)
	want := "@read/1\nAC\n+\n\"\"\n@read/1\nAC\n+\n\"\"\n@read/1\nAC\n+\n\"\"\n"
	if string(got) != want {
		t.Errorf("records = %q, want %q", string(got), want)
	}
}

func TestFastqWriter_EmptyRecordSet(t *testing.T) {
	path := filepath.Join(t.TempDir(), "empty.fq.gz")
	fw, err := NewFastqWriter(path, "/1", 1, 1)
	if err != nil {
		t.Fatalf("NewFastqWriter: %v", err)
	}
	if err := fw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	got := readGzipFile(t, path)
	if len(got) != 0 {
		t.Errorf("expected empty output, got %q", string(got))
	}
}

func TestNewFastqWriter_InvalidPath(t *testing.T) {
	_, err := NewFastqWriter(filepath.Join(t.TempDir(), "nonexistent-dir", "out.fq.gz"), "/1", 1, 1)
	if err == nil {
		t.Fatalf("expected an error creating a file in a nonexistent directory")
	}
}

func TestWriteInt(t *testing.T) {
	var buf bytes.Buffer
	w := bufio.NewWriter(&buf)
	tests := []struct {
		n    int
		want string
	}{
		{0, "0"},
		{7, "7"},
		{42, "42"},
		{1234567890, "1234567890"},
	}
	for _, tt := range tests {
		buf.Reset()
		if err := writeInt(w, tt.n); err != nil {
			t.Fatalf("writeInt(%d): %v", tt.n, err)
		}
		w.Flush()
		if got := buf.String(); got != tt.want {
			t.Errorf("writeInt(%d) = %q, want %q", tt.n, got, tt.want)
		}
	}
}
