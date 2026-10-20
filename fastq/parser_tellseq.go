package fastq

import (
	"github.com/shenwei356/bio/seqio/fastx"
)

// ---- tellseq -----------------------------

// findTellseqBarcode locates a ":<ACGTN+>" suffix in id (TELLseq's inline
// barcode format), optionally followed by "-" and a suffix (any
// non-whitespace content, e.g. Cell Ranger's "-1" GEM well convention).
// Candidate colons are tried right to left, mirroring the original
// right-anchored match (id, containing no whitespace, can only match at
// the very end) when there's no such suffix. Returns the barcode's
// [start,end) -- including any such suffix -- and the index of the
// leading ':' to splice from.
func findTellseqBarcode(id []byte) (start, end, colon int, ok bool) {
	n := len(id)
	for i := n - 1; i >= 0; i-- {
		if id[i] != ':' {
			continue
		}
		tail := id[i+1:]
		j := 0
		for j < len(tail) && isBase(tail[j]) {
			j++
		}
		if j == 0 {
			continue
		}
		rest := tail[j:]
		if len(rest) == 0 || (rest[0] == '-' && len(rest) > 1) {
			return i + 1, n, i, true
		}
	}
	return 0, 0, 0, false
}

// Tellseq2Corefq parses a TELLseq-format record — barcode inline at the end
// of the ID as ":<bases>" — into a CoreFq. ok is false when no such suffix
// is present
func Tellseq2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	start, end, colon, found := findTellseqBarcode(id)
	if !found {
		return CoreFq{
			ID:       id,
			Seq:      rec.Seq.Seq,
			Qual:     rec.Seq.Qual,
			CASAVA:   casava,
			Barcode:  nil,
			Valid:    false,
			Comments: desc,
		}, false
	}
	bc := append([]byte(nil), id[start:end]...)
	id = append(id[:colon], id[end:]...)

	return CoreFq{
		ID:       id,
		Seq:      rec.Seq.Seq,
		Qual:     rec.Seq.Qual,
		CASAVA:   casava,
		Barcode:  bc,
		Valid:    isValidNuc(bc),
		Comments: desc,
	}, true
}
