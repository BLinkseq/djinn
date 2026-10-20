package fastq

import (
	"bufio"
	"bytes"
	"fmt"
	"strconv"

	"github.com/biogo/hts/sam"
)

// FASTQ text format primitives used by the ToXxx writers below.
const FwSlash = '/'
const TabSep = '\t'
const Newline = '\n'
const PlusSign = '+'
const FastqAt = '@'
const Colon = ':'

var BxTagConst = []byte{'B', 'X', ':', 'Z', ':'}
var VxTagConst = []byte{'V', 'X', ':', 'i', ':'}

// CoreFq is a zero-copy, per-record view over the parts of a FASTQ record
// that matter for linked-read barcode handling: sequence/quality, the
// paired-end (CASAVA) marker, the isolated barcode, and (for direct BAM
// output) whatever SAM-style aux text is left in the description once the
// barcode and CASAVA marker have been pulled out of it.
//
// Every []byte field aliases memory owned by the *fastx.Record (or, more
// precisely, the fastx.Reader that produced it — see below) CoreFq was
// built from; building one does not allocate beyond copying out the short
// Barcode/CASAVA substrings that must survive being spliced out of ID/Desc.
//
// Lifetime: fastx.Reader.Read() reuses a single Record and its internal
// buffers across calls — it does not allocate a fresh one each time. A
// CoreFq is therefore only valid until the next call to that reader's
// Read(): build it, use it (convert/write it), and let it go before
// reading the next record. Never hand a CoreFq, or any of its fields, to
// another goroutine while the reader keeps advancing — that reintroduces
// the same class of data corruption bug this design is built to avoid (see
// HaplotagBX's bc-aliasing fix in processfq.go). For real concurrency, fan
// out over fastx.Reader.ChunkChan(bufferSize, chunkSize) instead: it calls
// Record.Clone() before a record crosses the channel, so each chunk owns
// independent memory and a CoreFq built from it inside a worker goroutine
// is safe on its own.
type CoreFq struct {
	ID       []byte // read ID, CASAVA suffix removed
	Seq      []byte // nucleotide sequence
	Qual     []byte // PHRED quality string
	CASAVA   []byte // paired-end marker: "1"/"2" (old style) or "1:N:0:ATCG" (new style); nil if none was found
	Barcode  []byte // isolated linked-read barcode; nil if none was found
	Valid    bool   // whether Barcode passes the Invalid pattern; meaningless when Barcode is nil
	Comments []byte // remaining description text (SAM-style aux tokens) once CASAVA/Barcode are removed
}

var CasavaConst = []byte{':', 'N', ':', '1', '8', ':', 'A', 'T', 'C', 'G'}

// ----- shared scanning helpers -------------------------------------------
//
// These replace regexp.FindSubmatchIndex-based matching (IlluminaNew,
// IlluminaOld, StdBx, Tellseq, Stlfr in barcodes.go) with single-pass byte
// scans that match the exact same rules. Regexp is fine at the call volumes
// those functions see elsewhere, but CoreFq is meant to be built once per
// record at tens-of-millions-of-records scale, so it avoids the engine
// overhead entirely.

func isSpace(b byte) bool {
	return b == ' ' || b == '\t'
}

func isBase(b byte) bool {
	switch b {
	case 'A', 'C', 'G', 'T', 'N':
		return true
	}
	return false
}

func isDigit(b byte) bool {
	return b >= '0' && b <= '9'
}

func isAlphaNumeric(b byte) bool {
	return b >= 'A' && b <= 'Z' || b >= 'a' && b <= 'z' || isDigit(b)
}

// nextField returns the [start,end) range of the next whitespace-delimited
// field in data at or after idx, and the offset to resume scanning from.
// ok is false once no more fields remain.
func nextField(data []byte, idx int) (start, end, next int, ok bool) {
	n := len(data)
	for idx < n && isSpace(data[idx]) {
		idx++
	}
	if idx >= n {
		return 0, 0, n, false
	}
	start = idx
	for idx < n && !isSpace(data[idx]) {
		idx++
	}
	return start, idx, idx, true
}

// findField scans desc field-by-field for one that starts with prefix, and
// returns the [valStart,valEnd) range of the value after prefix along with
// the [fieldStart,fieldEnd) range of the whole field (for splicing it back
// out of desc). ok is false if no field has that prefix.
func findField(desc, prefix []byte) (valStart, valEnd, fieldStart, fieldEnd int, ok bool) {
	for i := 0; ; {
		fs, fe, next, more := nextField(desc, i)
		if !more {
			return 0, 0, 0, 0, false
		}
		i = next
		if fe-fs > len(prefix) && bytes.HasPrefix(desc[fs:fe], prefix) {
			return fs + len(prefix), fe, fs, fe, true
		}
	}
}

// spliceField removes data[start:end) in place, absorbing one adjacent
// separator byte if present so a stray double-separator isn't left behind.
func spliceField(data []byte, start, end int) []byte {
	if end < len(data) && isSpace(data[end]) {
		end++
	} else if start > 0 && isSpace(data[start-1]) {
		start--
	}
	return append(data[:start], data[end:]...)
}

// isIlluminaNewField reports whether f matches the new-style CASAVA marker
// "[12]:[YN]:<digits>:<alnum+>" (e.g. "1:N:0:ATCG"), mirroring IlluminaNew.
func isIlluminaNewField(f []byte) bool {
	if len(f) < 7 || (f[0] != '1' && f[0] != '2') || f[1] != ':' {
		return false
	}
	if f[2] != 'Y' && f[2] != 'N' {
		return false
	}
	if f[3] != ':' {
		return false
	}
	i := 4
	digitsStart := i
	for i < len(f) && isDigit(f[i]) {
		i++
	}
	if i == digitsStart || i >= len(f) || f[i] != ':' {
		return false
	}
	i++
	if i >= len(f) {
		return false
	}
	for ; i < len(f); i++ {
		if !isAlphaNumeric(f[i]) {
			return false
		}
	}
	return true
}

// findIlluminaNew locates a new-style CASAVA field in desc, mirroring
// IlluminaNew.FindIndex(desc).
func findIlluminaNew(desc []byte) (start, end int, ok bool) {
	for i := 0; ; {
		fs, fe, next, more := nextField(desc, i)
		if !more {
			return 0, 0, false
		}
		i = next
		if isIlluminaNewField(desc[fs:fe]) {
			return fs, fe, true
		}
	}
}

// findIlluminaOld locates an old-style "/1" or "/2" CASAVA suffix at the
// end of id, mirroring IlluminaOld.FindIndex(id) — id never contains
// whitespace (fastx.Reader splits it from desc at the first whitespace),
// so the only place IlluminaOld's "/[12](?:\s|$)" can match is the very end.
func findIlluminaOld(id []byte) (start, end int, ok bool) {
	n := len(id)
	if n < 2 || id[n-2] != '/' {
		return 0, 0, false
	}
	if id[n-1] != '1' && id[n-1] != '2' {
		return 0, 0, false
	}
	return n - 2, n, true
}

// splitCASAVA removes the paired-end/CASAVA markers from id and desc -- new style
// ("1:N:0:ATCG") lives in desc, old style ("/1" or "/2") at the end of id -- and returns
// the shortened id and desc plus a copy of the marker that was removed.
//
// Both are always taken out, so a read carrying both (e.g. "read#1_2_3/1 1:N:0:ATCG") still
// leaves a clean id for the barcode search that follows. When both are present the new-style
// marker is the one returned, since it says more.
func splitCASAVA(id, desc []byte) (newID, newDesc, casava []byte) {
	newID, newDesc = id, desc
	// The markers are copied out before any splice: the splices reuse the same backing arrays,
	// which would otherwise overwrite what casava points at.
	if s, e, ok := findIlluminaOld(newID); ok {
		casava = append([]byte(nil), newID[s+1:e]...) // digit only, "/" dropped
		newID = newID[:s]                             // the marker is the very end of the id
	}
	if s, e, ok := findIlluminaNew(newDesc); ok {
		casava = append([]byte(nil), newDesc[s:e]...)
		newDesc = spliceField(newDesc, s, e)
	}
	return newID, newDesc, casava
}

// ---- Conversions -------------------------------------------------

// Convert the fastq record to haplotagging format and write to the write buffer
func (rec *CoreFq) ToHaplotagging(writer *bufio.Writer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
	// CASAVA may be empty (no "/1" or "N:0:ATCG"-style marker was present
	// on the source record) -- omit the suffix rather than index
	// CASAVA[0] unconditionally, which would panic on exactly that input.
	// Same guard as ToStandard.
	if len(rec.CASAVA) > 0 {
		writer.WriteByte(FwSlash)
		writer.WriteByte(rec.CASAVA[0])
	}
	writer.WriteByte(TabSep)
	//TODO ONLY WRITE SAM-COMPLIANT COMMENTS
	writer.Write(BxTagConst)
	writer.Write(rec.Barcode)
	writer.WriteByte(Newline)
	writer.Write(rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to standard format and write to the write buffer.
// This is the generic passthrough writer used by callers that only
// filter/classify records (count, extract, invalid, sample, singletons)
// rather than converting between linked-read technologies: output is always
// the same BX:Z:/VX:i: tag form regardless of the record's original source
// format.
func (rec *CoreFq) ToStandard(writer *bufio.Writer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
	// CASAVA may be empty (no "/1" or "N:0:ATCG"-style marker was present
	// on the source record) -- omit the suffix rather than index CASAVA[0]
	// unconditionally, which would panic on exactly that input.
	if len(rec.CASAVA) > 0 {
		writer.WriteByte(FwSlash)
		writer.WriteByte(rec.CASAVA[0])
	}
	writer.WriteByte(TabSep)
	writer.Write(BxTagConst)
	writer.Write(rec.Barcode)
	writer.WriteByte(TabSep)
	writer.Write(VxTagConst)
	if rec.Valid {
		writer.WriteByte('1')
	} else {
		writer.WriteByte('0')
	}
	// Comments holds whatever other SAM-style aux text remained in the
	// description once CASAVA/the barcode were pulled out of it -- carry
	// it through rather than silently dropping it.
	if len(rec.Comments) > 0 {
		writer.WriteByte(TabSep)
		writer.Write(rec.Comments)
	}
	writer.WriteByte(Newline)
	writer.Write(rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to tellseq format and write to the write buffer
func (rec *CoreFq) ToTellseq(writer *bufio.Writer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
	writer.WriteByte(Colon)
	writer.Write(rec.Barcode)
	writer.WriteByte(TabSep)
	if len(rec.CASAVA) == 1 {
		// then the CASAVA is a single byte, like from /1 or /2
		// only use the first byte
		writer.WriteByte(rec.CASAVA[0])
		writer.Write(CasavaConst)
	} else {
		writer.Write(rec.CASAVA)
	}
	writer.WriteByte(Newline)
	writer.Write(rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to haplotagging format and write it the write buffer
func (rec *CoreFq) ToStlfr(writer *bufio.Writer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
	writer.WriteByte('#')
	writer.Write(rec.Barcode)
	writer.WriteByte(TabSep)
	if n := len(rec.CASAVA); n == 1 {
		// then the CASAVA is a single byte, like from /1 or /2
		// only use the first byte
		writer.WriteByte(rec.CASAVA[0])
		writer.Write(CasavaConst)
	} else {
		writer.Write(rec.CASAVA)
	}
	writer.WriteByte(Newline)
	writer.Write(rec.Seq)
	writer.WriteByte(Newline)
	writer.WriteByte(PlusSign)
	writer.WriteByte(Newline)
	writer.Write(rec.Qual)
	writer.WriteByte(Newline)
}

// Convert the fastq record to 10X format and write to the write buffer
func (rec *CoreFq) ToTenX(writer *bufio.Writer) {
	writer.WriteByte(FastqAt)
	writer.Write(rec.ID)
	writer.WriteByte(Colon)
	writer.Write(rec.Barcode)
	writer.WriteByte(TabSep)
	// only use the first byte of CASAVA, which will either be 1 or 2
	if n := len(rec.CASAVA); n == 1 {
		writer.WriteByte(rec.CASAVA[0])
		writer.Write(CasavaConst)
	} else {
		writer.Write(rec.CASAVA)
	}
	writer.WriteByte(Newline)
	// only write barcode as R1 seq prefix. CASAVA may be empty (no "/1"
	// or "N:0:ATCG"-style marker was present, e.g. single-end input) --
	// indexing CASAVA[0] unconditionally would panic on exactly that
	// input, same bug ToHaplotagging/ToStandard already guard against.
	// Default to the R1 (barcode-bearing) branch when the mate can't be
	// told apart: 10X's barcode only ever lives on the "first" read, so
	// a lone single-end record is that read by definition.
	if len(rec.CASAVA) == 0 || rec.CASAVA[0] == '1' {
		writer.Write(rec.Barcode)
		writer.Write(rec.Seq)
		writer.WriteByte(Newline)
		writer.WriteByte(PlusSign)
		writer.WriteByte(Newline)
		for range rec.Barcode {
			writer.WriteByte('I')
		}
		writer.Write(rec.Qual)
	} else {
		writer.Write(rec.Seq)
		writer.WriteByte(Newline)
		writer.WriteByte(PlusSign)
		writer.WriteByte(Newline)
		writer.Write(rec.Qual)
	}
	writer.WriteByte(Newline)
}

// --- build a sam.Record for direct BAM output ----------------------------
//
// Writing via sam.Record + bam.Writer (rather than hand-rolling the BAM
// binary layout here) lets biogo/hts own the parts that are easy to get
// subtly wrong: the file-level BAM header/magic, 2-bit sequence packing,
// the index bin, and binary aux-field encoding.

// parseAuxFields converts SAM-text aux fields ("XX:i:123\tYY:Z:foo ...",
// tab or space separated) into sam.Aux values. Numeric types collapse to
// 'i' (int32) and 'f' (float32) for simplicity -- not size-optimal like
// htslib, but correct. B/H array tags are unsupported and return an error.
func parseAuxFields(comments []byte) ([]sam.Aux, error) {
	if len(comments) == 0 {
		return nil, nil
	}
	fields := bytes.FieldsFunc(comments, func(r rune) bool { return r == '\t' || r == ' ' })
	aux := make([]sam.Aux, 0, len(fields))
	for _, f := range fields {
		if len(f) < 5 || f[2] != ':' || f[4] != ':' {
			return nil, fmt.Errorf("bad aux field: %q", f)
		}
		tag, typ, val := sam.NewTag(string(f[0:2])), f[3], f[5:]
		var a sam.Aux
		var err error
		switch typ {
		case 'A':
			a, err = sam.NewAux(tag, sam.ASCII(val[0]))
		case 'i':
			var n int64
			n, err = strconv.ParseInt(string(val), 10, 32)
			if err == nil {
				a, err = sam.NewAux(tag, int32(n))
			}
		case 'f':
			var fl float64
			fl, err = strconv.ParseFloat(string(val), 32)
			if err == nil {
				a, err = sam.NewAux(tag, float32(fl))
			}
		case 'Z':
			a, err = sam.NewAux(tag, sam.Text(val))
		case 'H':
			return nil, fmt.Errorf("H aux type unsupported: %q", tag)
		case 'B':
			return nil, fmt.Errorf("B aux type unsupported: %q", tag)
		default:
			return nil, fmt.Errorf("unknown aux type %c in %q", typ, tag)
		}
		if err != nil {
			return nil, err
		}
		aux = append(aux, a)
	}
	return aux, nil
}

var bxTag = sam.NewTag("BX")
var vxTag = sam.NewTag("VX")

// ToSamRecord builds an unmapped sam.Record from c, suitable for writing
// with a bam.Writer. All reads are unmapped: Ref/MateRef are nil, Pos/
// MatePos are -1, and both the Unmapped and MateUnmapped flags are always
// set (there is no mapped mate to report either way, paired or not).
//
// The barcode and its validity are written out as BX:Z:/VX:i: aux tags
// from c.Barcode/c.Valid, the same way ToStandard does for FASTQ output --
// the barcode parsers splice BX:Z: (and VX:i:) out of the description into
// those fields, so c.Comments alone never carries it, and a writer that
// only serialized Comments (as this one used to) would silently drop the
// barcode from the BAM output entirely.
func (c *CoreFq) ToSamRecord() (*sam.Record, error) {
	rest, err := parseAuxFields(c.Comments)
	if err != nil {
		return nil, err
	}
	bx, err := sam.NewAux(bxTag, sam.Text(c.Barcode))
	if err != nil {
		return nil, err
	}
	validFlag := int32(0)
	if c.Valid {
		validFlag = 1
	}
	vx, err := sam.NewAux(vxTag, validFlag)
	if err != nil {
		return nil, err
	}
	aux := append([]sam.Aux{bx, vx}, rest...)

	qual := c.Qual
	if len(qual) == len(c.Seq) {
		converted := make([]byte, len(qual))
		for i, q := range qual {
			converted[i] = q - 33
		}
		qual = converted
	} else {
		qual = nil
	}

	rec, err := sam.NewRecord(string(c.ID), nil, nil, -1, -1, 0, 0, nil, c.Seq, qual, aux)
	if err != nil {
		return nil, err
	}

	rec.Flags = sam.Unmapped | sam.MateUnmapped
	if len(c.CASAVA) > 0 {
		rec.Flags |= sam.Paired
		switch c.CASAVA[0] {
		case '1':
			rec.Flags |= sam.Read1
		case '2':
			rec.Flags |= sam.Read2
		}
	}
	return rec, nil
}
