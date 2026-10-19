package fastq

import "github.com/shenwei356/bio/seqio/fastx"

// ----- haplotagging --------------------------------------------------
// Haplotag2Corefq parses a haplotagging-format record — barcode inline as a
// "BX:Z:<bc>" tag in the description — into a CoreFq. ok is false when no
// BX:Z: tag is present.
//
// Validity is inferred from the barcode text (see Invalid) unless an
// explicit, well-formed "VX:i:0" or "VX:i:1" field is also present, in
// which case that value takes precedence -- same honor-VX-if-present
// precedence as xam.FindBarcode.
func Haplotag2Corefq(rec *fastx.Record) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)

	valStart, valEnd, fs, fe, found := findField(desc, BXTAG)
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
	bc := append([]byte(nil), desc[valStart:valEnd]...)
	desc = spliceField(desc, fs, fe)

	valid := !isInvalidBarcode(bc)
	if vs, ve, vfs, vfe, vfound := findField(desc, VXTAG); vfound && ve-vs == 1 && (desc[vs] == '0' || desc[vs] == '1') {
		valid = desc[vs] == '1'
		desc = spliceField(desc, vfs, vfe)
	}

	return CoreFq{
		ID:       id,
		Seq:      rec.Seq.Seq,
		Qual:     rec.Seq.Qual,
		CASAVA:   casava,
		Barcode:  bc,
		Valid:    valid,
		Comments: desc,
	}, true
}
