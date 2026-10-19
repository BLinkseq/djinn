package xam

import (
	"encoding/binary"
	"strings"

	"github.com/biogo/hts/sam"
)

// VX sam tag
var VxTag = sam.Tag{'V', 'X'}

// BX sam tag
var BxTag = sam.Tag{'B', 'X'}

// Returns a true if a barcode is valid (i.e. not invalid) in either
// haplotagging, tellseq, or stlfr formats. Only the portion before an
// optional "-"+suffix (see hasHyphenSuffix in parsers.go) is checked --
// no recognized barcode shape contains a hyphen itself, so any hyphen
// present marks the start of a technology-agnostic trailing tag (e.g. a
// GEM well index) that is never itself part of the barcode's validity.
func IsValid(barcode string) bool {
	if i := strings.IndexByte(barcode, '-'); i != -1 {
		barcode = barcode[:i]
	}
	if strings.IndexByte(barcode, 'N') != -1 {
		return false
	}
	if strings.HasPrefix(barcode, "0_") || strings.HasSuffix(barcode, "_0") || strings.Contains(barcode, "_0_") {
		return false
	}
	for _, p := range [4]string{"A00", "B00", "C00", "D00"} {
		if strings.Contains(barcode, p) {
			return false
		}
	}
	return true
}

// isBase reports whether b is one of the DNA bases used in barcode text
// (including the ambiguity code N).
func isBase(b byte) bool {
	switch b {
	case 'A', 'C', 'G', 'T', 'N':
		return true
	}
	return false
}

func isDigit(b byte) bool { return b >= '0' && b <= '9' }

// isStlfrBody reports whether s matches "[0-9]+_[0-9]+_[0-9]+", optionally
// followed by "-" and a suffix (see hasHyphenSuffix in parsers.go).
func isStlfrBody(s string) bool {
	i, parts := 0, 0
	for parts < 3 {
		start := i
		for i < len(s) && isDigit(s[i]) {
			i++
		}
		if i == start {
			return false
		}
		parts++
		if parts == 3 {
			break
		}
		if i >= len(s) || s[i] != '_' {
			return false
		}
		i++
	}
	return hasHyphenSuffix(s[i:])
}

// findInlineBarcode locates a linked-read barcode inline in name, in
// either TELLseq (":<bases>" suffix) or stLFR ("#<n>_<n>_<n>" suffix)
// style, mirroring what the pattern `(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`
// matched, but as a pair of native byte scans instead of a regexp. Either
// shape may be followed by "-" and a suffix (see hasHyphenSuffix); the
// returned barcode includes that suffix, same as a BX:Z tag's value
// would never have had it stripped in the first place.
func findInlineBarcode(name string) (barcode string, ok bool) {
	// TELLseq: "<bases>" (optionally "-"+suffix) right after a ':',
	// reaching to the end of name. Scan candidate colons from the right,
	// mirroring the original right-anchored match when there's no suffix.
	for i := len(name) - 1; i >= 0; i-- {
		if name[i] != ':' {
			continue
		}
		tail := name[i+1:]
		j := 0
		for j < len(tail) && isBase(tail[j]) {
			j++
		}
		if j > 0 && hasHyphenSuffix(tail[j:]) {
			return tail, true
		}
	}

	// stLFR: "#<n>_<n>_<n>" (optionally "-"+suffix) at the end of name.
	// Scan every '#' left to right (mirroring the leftmost-match behavior
	// of the regexp) and take the first whose remainder matches.
	for j := range name {
		if name[j] != '#' {
			continue
		}
		if body := name[j+1:]; isStlfrBody(body) {
			return body, true
		}
	}
	return "", false
}

// Search for an return a linked read barcode, ("" if not found) and the value of the VX tag (bool).
// First searches for a BX tag, and if that isn't found, searches the record ID for a tellseq/stlfr style barcode.
// If nothing was found, returns ("", false). If a barcode was identified and
// a VX tag wasnt, the VX will be inferred from the barcode.
func FindBarcode(rec *sam.Record) (string, bool) {
	bxVal, hasBX := GetStringTag(rec, "BX")
	vxVal, hasVX := GetVX(rec)
	if !hasBX {
		bxVal, hasBX = findInlineBarcode(rec.Name)
	}
	if !hasVX && bxVal != "" {
		vxVal = IsValid(bxVal)
	}
	return bxVal, vxVal
}

// SetBX sets a string aux tag on a record
func SetBxString(rec *sam.Record, val *string) {
	aux := make(sam.Aux, 3+len(*val))
	aux[0], aux[1], aux[2] = 'B', 'X', 'Z'
	copy(aux[3:], *val)

	for i, a := range rec.AuxFields {
		if a.Tag() == BxTag {
			rec.AuxFields[i] = aux
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, aux)
}

// SetBX sets a string aux tag on a record
func SetBxByte(rec *sam.Record, val *[]byte) {
	aux := make(sam.Aux, 3+len(*val))
	aux[0], aux[1], aux[2] = 'B', 'X', 'Z'
	copy(aux[3:], *val)

	for i, a := range rec.AuxFields {
		if a.Tag() == BxTag {
			rec.AuxFields[i] = aux
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, aux)
}

/*
Convenience function to return the VX tag as a bool.
Returns (bool, bool), where the 2nd bool is whether
the VX tag was present. i.e.:

true, true -> VX:i = 1 and it was present

true, false -> VX:i = 1 and it wasn't present (impossible)

false, true -> VX:i = 0 and it was present

false, false VX:i = 0 because it wasn't present
*/
func GetVX(rec *sam.Record) (bool, bool) {
	vx, hasVX := GetIntTag(rec, "VX")
	vxVal := vx == 1
	return vxVal, hasVX
}

// SetVX sets an integer (0/1) auxiliary tag on a record.
func SetVX(rec *sam.Record, isValid bool) {
	var v int32
	if isValid {
		v = 1
	}
	buf := make([]byte, 4)
	binary.LittleEndian.PutUint32(buf, uint32(v))

	aux := sam.Aux{'V', 'X', 'i', 0, 0, 0, 0}
	copy(aux[3:7], buf)

	for i, a := range rec.AuxFields {
		if a.Tag() == VxTag {
			// Replace the whole field rather than copying into the
			// existing byte range: an aux field parsed from SAM text
			// (e.g. "VX:i:1") may have been encoded by the library as
			// a narrower type ('c'/'C'/'s'/'S', 1-2 bytes) instead of
			// the 4-byte 'i' this function assumes, and copy(a[3:7], buf)
			// would panic with an out-of-range slice on such a field.
			rec.AuxFields[i] = aux
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, aux)
}

// Return the value of the string tag `tag` (XX:Z) for a sam.Record
func GetStringTag(r *sam.Record, tag string) (string, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			if s, ok := aux.Value().(string); ok {
				return s, true
			}
		}
	}
	return "", false
}

// Return the value of the integer tag `tag` (XX:i) for a sam.Record
func GetIntTag(r *sam.Record, tag string) (int, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			// sam.Aux.Value() never returns a plain "int" for integer-kind
			// tags: depending on the field's stored width it returns
			// int8/uint8/int16/uint16/int32/uint32, so a bare .(int)
			// assertion never succeeds and this always reported "not found".
			switch v := aux.Value().(type) {
			case int8:
				return int(v), true
			case uint8:
				return int(v), true
			case int16:
				return int(v), true
			case uint16:
				return int(v), true
			case int32:
				return int(v), true
			case uint32:
				return int(v), true
			}
		}
	}
	return 0, false
}

// AlignedEnd is the position just past the last reference base a read covers, so Pos and AlignedEnd
// give the stretch of the genome it spans. Only CIGAR operations that walk along the reference count
// (matches, deletions, skips); a CIGAR that covers no reference at all counts as one base, as htslib's
// bam_endpos does.
func AlignedEnd(r *sam.Record) int {
	n := 0
	for _, op := range r.Cigar {
		switch op.Type() {
		case sam.CigarMatch, sam.CigarDeletion, sam.CigarSkipped, sam.CigarEqual, sam.CigarMismatch:
			n += op.Len()
		}
	}
	if n == 0 {
		n = 1
	}
	return r.Pos + n
}
