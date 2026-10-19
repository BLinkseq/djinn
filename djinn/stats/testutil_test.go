package stats

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

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
	return <-outCh, fnErr
}

// captureStderr is captureStdout for os.Stderr.
func captureStderr(t *testing.T, fn func() error) (string, error) {
	t.Helper()
	orig := os.Stderr
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	os.Stderr = w
	outCh := make(chan string, 1)
	go func() {
		var buf bytes.Buffer
		io.Copy(&buf, r)
		outCh <- buf.String()
	}()

	fnErr := fn()

	w.Close()
	os.Stderr = orig
	return <-outCh, fnErr
}

func tempPath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join(t.TempDir(), name)
}

const samHeader = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000000\n@SQ\tSN:c2\tLN:1000000\n"

// samLine is one unpaired 10M read; tags are tab-joined, e.g. "BX:Z:AAAA\tVX:i:1".
func samLine(name, contig string, pos int, tags string) string {
	return samLineCigar(name, 0, contig, pos, "10M", 10, tags)
}

func samLineCigar(name string, flag int, contig string, pos int, cigar string, seqLen int, tags string) string {
	line := fmt.Sprintf("%s\t%d\t%s\t%d\t60\t%s\t*\t0\t0\t%s\t%s", name, flag, contig, pos, cigar,
		strings.Repeat("A", seqLen), strings.Repeat("I", seqLen))
	if tags != "" {
		line += "\t" + tags
	}
	return line + "\n"
}

// writeSam writes a SAM file: the standard header followed by the given records.
func writeSam(t *testing.T, records ...string) string {
	t.Helper()
	path := tempPath(t, "in.sam")
	if err := os.WriteFile(path, []byte(samHeader+strings.Join(records, "")), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// statsOf runs StatsXam on path and returns the output split into lines (header excluded).
func statsOf(t *testing.T, path string) []string {
	t.Helper()
	out, err := captureStdout(t, func() error { return StatsXam(path, 1) })
	if err != nil {
		t.Fatalf("StatsXam: %v", err)
	}
	lines := strings.Split(strings.TrimSuffix(out, "\n"), "\n")
	if lines[0] != strings.TrimSuffix(header, "\n") {
		t.Fatalf("missing header, got %q", lines[0])
	}
	return lines[1:]
}

func sortedLines(s string) []string {
	lines := strings.Split(strings.TrimSuffix(s, "\n"), "\n")
	sort.Strings(lines)
	return lines
}
