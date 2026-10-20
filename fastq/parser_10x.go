package fastq

import (
	"bytes"

	"github.com/shenwei356/bio/seqio/fastx"
)

func isValidNuc(seq []byte) bool {
	if bytes.IndexByte(seq, 'N') != -1 {
		return false // Found 'N'
	}
	return true
}

// Parse an R1 10X style FASTQ record to identify the barcode against a 10X barcode set (bc_map), returning
// a CoreFq object where the barcode, if found, is removed removed and put into the Barcode field of the struct.
func TenXR12Corefq(rec *fastx.Record, bc_map map[string]struct{}) (core CoreFq, ok bool) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)
	bc := rec.Seq.Seq[0:16]
	if _, ok := bc_map[string(bc)]; ok {
		bc = append([]byte(nil), rec.Seq.Seq[0:16]...)
		rec.Seq.Seq = append([]byte(nil), rec.Seq.Seq[16:]...)
		rec.Seq.Qual = append([]byte(nil), rec.Seq.Qual[16:]...)
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

// Parse an R2 10X style FASTQ record into a CoreFq struct. Since R2 reads don't carry barcodes, only the CASAVA
// is split from the read.
func TenXR22Corefq(rec *fastx.Record) (core CoreFq) {
	id, desc, casava := splitCASAVA(rec.ID, rec.Desc)
	return CoreFq{
		ID:       id,
		Seq:      rec.Seq.Seq,
		Qual:     rec.Seq.Qual,
		CASAVA:   casava,
		Barcode:  nil,
		Valid:    false,
		Comments: desc,
	}
}
