package barcodes

import "iter"

var invalidHaplotag = [12]byte{'A', '0', '0', 'C', '0', '0', 'B', '0', '0', 'D', '0', '0'}

func writeHaplotagPart(buf *[12]byte, off int, prefix byte, n int) {
	buf[off] = prefix
	buf[off+1] = byte('0' + n/10)
	buf[off+2] = byte('0' + n%10)
}

func haplotagSeq() iter.Seq[[12]byte] {
	return func(yield func([12]byte) bool) {
		var buf [12]byte
		for a := 1; a <= 96; a++ {
			writeHaplotagPart(&buf, 0, 'A', a)
			for c := 1; c <= 96; c++ {
				writeHaplotagPart(&buf, 3, 'C', c)
				for b := 1; b <= 96; b++ {
					writeHaplotagPart(&buf, 6, 'B', b)
					for d := 1; d <= 96; d++ {
						writeHaplotagPart(&buf, 9, 'D', d)
						if !yield(buf) {
							return
						}
					}
				}
			}
		}
	}
}

type Haplotagging struct {
	invalid [12]byte
	next    func() ([12]byte, bool)
	stop    func()
}

func NewHaplotagging() *Haplotagging {
	next, stop := iter.Pull(haplotagSeq())
	return &Haplotagging{invalid: invalidHaplotag, next: next, stop: stop}
}

func (h *Haplotagging) Next() ([12]byte, bool) { return h.next() }

// func (h *Haplotagging) Invalid() [12]byte      { return h.invalid }
func (h *Haplotagging) Close() { h.stop() }

// func (h *Haplotagging) InvalidInto(dst []byte) int { return copy(dst, h.invalid[:]) }
func (h *Haplotagging) GetInvalid() []byte { return h.invalid[:12] }
func (h *Haplotagging) MaxLen() int        { return 12 }
func (h *Haplotagging) NextInto(dst []byte) (int, bool) {
	b, ok := h.next()
	if !ok {
		return 0, false
	}
	return copy(dst, b[:]), true
}
