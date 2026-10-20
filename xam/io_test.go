package xam

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/biogo/hts/sam"
)

func newHeader(t *testing.T) *sam.Header {
	t.Helper()
	h, err := sam.NewHeader(nil, nil)
	if err != nil {
		t.Fatalf("sam.NewHeader: %v", err)
	}
	return h
}

func TestNewPGFirstProgram(t *testing.T) {
	h := newHeader(t)
	pg := NewPG(h, "djinn sam convert haplotagging in.bam")
	if pg.UID() != "djinn" {
		t.Errorf("UID = %q, want %q", pg.UID(), "djinn")
	}
	if pg.Previous() != "" {
		t.Errorf("Previous = %q, want empty (no prior programs)", pg.Previous())
	}
}

func TestNewPGChainsToTail(t *testing.T) {
	h := newHeader(t)
	if err := h.AddProgram(sam.NewProgram("bwa", "bwa", "bwa mem", "", "0.7.17")); err != nil {
		t.Fatalf("AddProgram: %v", err)
	}
	pg := NewPG(h, "djinn sam convert haplotagging in.bam")
	if pg.Previous() != "bwa" {
		t.Errorf("Previous = %q, want %q", pg.Previous(), "bwa")
	}
	if pg.UID() != "djinn" {
		t.Errorf("UID = %q, want %q", pg.UID(), "djinn")
	}
}

func TestNewPGUniquifiesUID(t *testing.T) {
	h := newHeader(t)
	if err := h.AddProgram(sam.NewProgram("djinn", "djinn", "djinn cmd1", "", "3.0")); err != nil {
		t.Fatalf("AddProgram: %v", err)
	}
	pg := NewPG(h, "djinn sam convert haplotagging in.bam")
	if pg.UID() != "djinn.1" {
		t.Errorf("UID = %q, want %q", pg.UID(), "djinn.1")
	}

	if err := h.AddProgram(pg); err != nil {
		t.Fatalf("AddProgram: %v", err)
	}
	pg2 := NewPG(h, "djinn sam convert haplotagging in.bam")
	if pg2.UID() != "djinn.2" {
		t.Errorf("UID = %q, want %q", pg2.UID(), "djinn.2")
	}
}

func TestFileOrStdinWithFile(t *testing.T) {
	got, err := FileOrStdin("some.bam")
	if err != nil {
		t.Fatalf("FileOrStdin returned error: %v", err)
	}
	if got != "some.bam" {
		t.Errorf("FileOrStdin = %q, want %q", got, "some.bam")
	}
}

func TestBamNotStdoutSam(t *testing.T) {
	// When writing SAM (isSam == true) there should never be an error,
	// even when stdout is not redirected.
	if err := BamNotStdout(true); err != nil {
		t.Errorf("BamNotStdout(true) = %v, want nil", err)
	}
}

func TestXamReaderWriterRoundTripSAM(t *testing.T) {
	dir := t.TempDir()
	inPath := filepath.Join(dir, "in.sam")
	outPath := filepath.Join(dir, "out.sam")

	samText := "@HD\tVN:1.6\tSO:unsorted\n" +
		"@SQ\tSN:chr1\tLN:1000\n" +
		"read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
	if err := os.WriteFile(inPath, []byte(samText), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	recChan, reader := NewXamReaderChan(inPath, 10, IoBuf, 1)
	defer reader.Close()

	var recs []*sam.Record
	for rec := range recChan {
		recs = append(recs, rec)
	}
	if len(recs) != 1 {
		t.Fatalf("got %d records, want 1", len(recs))
	}
	if recs[0].Name != "read1" {
		t.Errorf("record name = %q, want %q", recs[0].Name, "read1")
	}

	hdr := reader.Header()
	writeChan, done := NewXamWriterChan(outPath, hdr, 10, IoBuf, 1, true)
	for _, rec := range recs {
		writeChan <- rec
	}
	close(writeChan)
	<-done

	out, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	if len(out) == 0 {
		t.Error("expected non-empty SAM output")
	}
}

func TestXamReaderEmptyChannelOnEmptyInput(t *testing.T) {
	dir := t.TempDir()
	inPath := filepath.Join(dir, "in.sam")
	samText := "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"
	if err := os.WriteFile(inPath, []byte(samText), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	recChan, reader := NewXamReaderChan(inPath, 10, IoBuf, 1)
	defer reader.Close()

	count := 0
	for range recChan {
		count++
	}
	if count != 0 {
		t.Errorf("got %d records, want 0", count)
	}
}
