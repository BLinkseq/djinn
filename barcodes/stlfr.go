package barcodes

import (
	"iter"
	"strconv"
)

const stlfrMaxLen = 14 // "1537_1537_1537"

type Barcode struct {
	Data [stlfrMaxLen]byte
	Len  int8
}

func (b Barcode) Bytes() []byte  { return b.Data[:b.Len] }
func (b Barcode) String() string { return string(b.Data[:b.Len]) }

var invalidStlfr = Barcode{Data: [stlfrMaxLen]byte{'0', '_', '0', '_', '0'}, Len: 5}

func stlfrSeq() iter.Seq[Barcode] {
	return func(yield func(Barcode) bool) {
		for x := 1; x < 1538; x++ {
			for y := 1; y < 1538; y++ {
				for z := 1; z < 1538; z++ {
					var bc Barcode
					n := bc.Data[:0]
					n = strconv.AppendInt(n, int64(x), 10)
					n = append(n, '_')
					n = strconv.AppendInt(n, int64(y), 10)
					n = append(n, '_')
					n = strconv.AppendInt(n, int64(z), 10)
					bc.Len = int8(len(n))
					if !yield(bc) {
						return
					}
				}
			}
		}
	}
}

type Stlfr struct {
	invalid Barcode
	next    func() (Barcode, bool)
	stop    func()
}

func NewStlfr() *Stlfr {
	next, stop := iter.Pull(stlfrSeq())
	return &Stlfr{invalid: invalidStlfr, next: next, stop: stop}
}

func (s *Stlfr) Next() (Barcode, bool) { return s.next() }

// func (s *Stlfr) Invalid() Barcode      { return s.invalid }
func (s *Stlfr) Close() { s.stop() }

func (s *Stlfr) NextInto(dst []byte) (int, bool) {
	bc, ok := s.next()
	if !ok {
		return 0, false
	}
	return copy(dst, bc.Bytes()), true
}
func (s *Stlfr) InvalidInto(dst []byte) int { return copy(dst, s.invalid.Bytes()) }
func (s *Stlfr) GetInvalid() []byte         { return s.invalid.Data[:s.invalid.Len] }
func (s *Stlfr) MaxLen() int                { return stlfrMaxLen }

/*
bcs := barcodes.NewHaplotagging()
/ or /
bcs := barcodes.NewStlfr()

defer bcs.Close()
bc, ok := bcs.Next()
if !ok {
	log.Fatal("exceeded max stlfr barcodes")
}
tag := bc.Bytes() // e.g. []byte("1_1_1")
*/
