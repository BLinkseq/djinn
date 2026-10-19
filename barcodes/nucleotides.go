package barcodes

import (
	"fmt"
	"iter"
)

const nucMaxLen = 24 // cap for the backing array

var nucAlphabet = [4]byte{'A', 'T', 'C', 'G'}

type NucBarcode struct {
	Data [nucMaxLen]byte
	Len  int
}

func (b NucBarcode) Bytes() []byte  { return b.Data[:b.Len] }
func (b NucBarcode) String() string { return string(b.Data[:b.Len]) }

func nucSeq(n int) iter.Seq[NucBarcode] {
	return func(yield func(NucBarcode) bool) {
		var bc NucBarcode
		bc.Len = n
		var rec func(pos int) bool
		rec = func(pos int) bool {
			if pos == n {
				return yield(bc) // copies bc by value, safe to keep mutating after
			}
			for _, ch := range nucAlphabet {
				bc.Data[pos] = ch
				if !rec(pos + 1) {
					return false
				}
			}
			return true
		}
		rec(0)
	}
}

type Nucleotides struct {
	n       int
	invalid NucBarcode
	next    func() (NucBarcode, bool)
	stop    func()
}

// NewTellseq creates a generic 18bp nucleotide barcode generator mimicking 10X (4^18 total combinations)
func NewTellseq() *Nucleotides {
	var invalid NucBarcode
	for i := range 18 {
		invalid.Data[i] = 'N'
	}
	invalid.Len = 18

	next, stop := iter.Pull(nucSeq(18))
	return &Nucleotides{n: 18, invalid: invalid, next: next, stop: stop}
}

// NewTenX creates a generic 16bp nucleotide barcode generator mimicking 10X (4^16 total combinations).
func NewTenX() *Nucleotides {
	var invalid NucBarcode
	for i := range 16 {
		invalid.Data[i] = 'N'
	}
	invalid.Len = 16

	next, stop := iter.Pull(nucSeq(16))
	return &Nucleotides{n: 16, invalid: invalid, next: next, stop: stop}
}

// NewGeneric creates a nucleotide barcode generator producing barcodes of length n
// over the ATCG alphabet (4^n total combinations).
func NewGeneric(n int) (*Nucleotides, error) {
	if n <= 0 || n > nucMaxLen {
		return nil, fmt.Errorf("barcodes: barcode n must be between 1 and %d, got %d", nucMaxLen, n)
	}

	var invalid NucBarcode
	for i := range n {
		invalid.Data[i] = 'N'
	}
	invalid.Len = n

	next, stop := iter.Pull(nucSeq(n))
	return &Nucleotides{n: n, invalid: invalid, next: next, stop: stop}, nil
}

func (t *Nucleotides) Next() (NucBarcode, bool) { return t.next() }

// func (t *Nucleotides) Invalid() NucBarcode      { return t.invalid }
func (t *Nucleotides) Close()             { t.stop() }
func (t *Nucleotides) GetInvalid() []byte { return t.invalid.Data[:t.n] }

// func (t *Nucleotides) InvalidInto(dst []byte) int { return copy(dst, t.invalid.Bytes()) }
func (t *Nucleotides) MaxLen() int { return t.n }
func (t *Nucleotides) NextInto(dst []byte) (int, bool) {
	bc, ok := t.next()
	if !ok {
		return 0, false
	}
	return copy(dst, bc.Bytes()), true
}
