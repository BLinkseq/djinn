package fastq

import "github.com/shenwei356/bio/seqio/fastx"

// isStlfrBody reports whether b matches "[0-9]+_[0-9]+_[0-9]+", optionally
// followed by "-" and a suffix (any non-whitespace content) -- a
// technology-agnostic allowance for a trailing tag, e.g. Cell Ranger's
// "-1" GEM well suffix convention, generalized here to every linked-read
// technology this package recognizes.
func isStlfrBody(b []byte) bool {
	i, parts := 0, 0
	for parts < 3 {
		start := i
		for i < len(b) && isDigit(b[i]) {
			i++
		}
		if i == start {
			return false
		}
		parts++
		if parts == 3 {
			break
		}
		if i >= len(b) || b[i] != '_' {
			return false
		}
		i++
	}
	rest := b[i:]
	return len(rest) == 0 || (rest[0] == '-' && len(rest) > 1)
}

// findStlfrBarcode locates a "#<n>_<n>_<n>" suffix in id, mirroring
// Stlfr.FindSubmatchIndex applied to id.
func findStlfrBarcode(id []byte) (start, end, hash int, ok bool) {
	for i := range id {
		if id[i] != '#' {
			continue
		}
		if body := id[i+1:]; isStlfrBody(body) {
			return i + 1, len(id), i, true
		}
	}
	return 0, 0, 0, false
}

// Stlfr2Corefq parses an stLFR-format record — barcode inline at the end of
// the ID as "#n_n_n" — into a CoreFq. ok is false when no such suffix is
// present
func Stlfr2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	start, end, hash, found := findStlfrBarcode(id)
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
	id = append(id[:hash], id[end:]...)

	return CoreFq{
		ID:       id,
		Seq:      rec.Seq.Seq,
		Qual:     rec.Seq.Qual,
		CASAVA:   casava,
		Barcode:  bc,
		Valid:    !isInvalidBarcode(bc),
		Comments: desc,
	}, true
}
