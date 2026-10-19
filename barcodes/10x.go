package barcodes

import (
	"bufio"
	"bytes"
	_ "embed"
	"fmt"
	"os"

	"github.com/klauspost/compress/zstd"
)

// The real 10X Genomics barcode whitelist ("4M-with-alts-february-2016"):
// one 16bp barcode per line, zstd-compressed (12.4MB, vs 15.8MB gzip -9 on
// the same data), embedded directly into the binary. 4,792,320 barcodes,
// all exactly 16 bases.
//
//go:embed 10x.zst
var tenXZst []byte

const tenXLen = 16

var invalidTenX = func() (b [tenXLen]byte) {
	for i := range b {
		b[i] = 'N'
	}
	return
}()

// TenXList is a Generator that streams the real, curated 10X barcode
// whitelist line by line (from the embedded asset), instead of
// enumerating every possible 16-mer the way Nucleotides/NewTenX does. Use
// this when assigning *new* 10X barcodes during a conversion -- e.g. the
// "10x" case in convert.fq.go's output-format switch.
type TenXList struct {
	dec     *zstd.Decoder
	scanner *bufio.Scanner
}

// NewTenXList opens a fresh stream over the embedded 10X barcode whitelist.
// Call Close when done to release the decoder. Expected to be used like
// NewHaplotagging()/NewStlfr()/NewTenX()
func NewTenXList() (*TenXList, error) {
	dec, err := zstd.NewReader(bytes.NewReader(tenXZst))
	if err != nil {
		return nil, fmt.Errorf("barcodes: opening embedded 10X whitelist: %w", err)
	}
	return &TenXList{dec: dec, scanner: bufio.NewScanner(dec)}, nil
}

// NextInto copies the next barcode in the whitelist into dst and reports
// whether one was available. As with every other Generator, dst must be
// used/consumed before the next call: bufio.Scanner's Bytes() (what this
// copies from) is only valid until the next Scan() anyway.
func (t *TenXList) NextInto(dst []byte) (int, bool) {
	if !t.scanner.Scan() {
		return 0, false
	}
	return copy(dst, t.scanner.Bytes()), true
}

func (t *TenXList) GetInvalid() []byte { return invalidTenX[:] }
func (t *TenXList) MaxLen() int        { return tenXLen }
func (t *TenXList) Close()             { t.dec.Close() }

// LoadTenXSet reads the full embedded 10X barcode whitelist into a
// membership set, for recognizing/validating an inline 10X barcode parsed
// out of a read's sequence (see fastq.TenX2Corefq's bc_map parameter) --
// a different access pattern than TenXList's sequential NextInto: this one
// needs the whole whitelist resident as a set for O(1) lookups, rather
// than handing out barcodes one at a time. Expected to be used for reading IN
// 10X data and parsing the first 16 bases.
func LoadTenXSet() (map[string]struct{}, error) {
	dec, err := zstd.NewReader(bytes.NewReader(tenXZst))
	if err != nil {
		return nil, fmt.Errorf("barcodes: opening embedded 10X whitelist: %w", err)
	}
	defer dec.Close()

	set := make(map[string]struct{}, 5_000_000)
	scanner := bufio.NewScanner(dec)
	for scanner.Scan() {
		// map keys always copy on insert (Go strings are immutable and a
		// []byte->string conversion here copies), so aliasing
		// scanner.Bytes() past this call is a non-issue.
		set[scanner.Text()] = struct{}{}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("barcodes: reading embedded 10X whitelist: %w", err)
	}
	return set, nil
}

// Read a 10X style barcode whitelist file into a map. Expects one barcode
// per line, does not check for correctness.
func LoadTenXFromFile(whitelist string) (map[string]struct{}, error) {
	dec, err := os.Open(whitelist)
	if err != nil {
		return nil, fmt.Errorf("failed opening 10X whitelist: %w", err)
	}
	defer dec.Close()

	set := make(map[string]struct{}, 5_000_000)
	scanner := bufio.NewScanner(dec)
	for scanner.Scan() {
		// map keys always copy on insert (Go strings are immutable and a
		// []byte->string conversion here copies), so aliasing
		// scanner.Bytes() past this call is a non-issue.
		set[scanner.Text()] = struct{}{}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("failed reading 10X whitelist: %w", err)
	}
	return set, nil
}
