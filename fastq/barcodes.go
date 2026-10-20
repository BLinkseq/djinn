package fastq

import "bytes"

// isInvalidBarcode reports whether bc is a known-bad haplotagging or
// stLFR barcode: it contains the ambiguity code 'N', an all-zero
// haplotagging segment (A00/B00/C00/D00), or an all-zero/leading-zero
// stLFR component (a leading "0_", a trailing "_0", or an internal
// "_0_"). Mirrors xam.IsValid's rule (inverted), just for a []byte
// barcode instead of a string -- including stripping an optional
// "-"+suffix before checking: no recognized barcode shape contains a
// hyphen itself, so any hyphen present marks the start of a
// technology-agnostic trailing tag (e.g. a GEM well index) that is never
// itself part of the barcode's validity.
func isInvalidBarcode(bc []byte) bool {
	if i := bytes.IndexByte(bc, '-'); i != -1 {
		bc = bc[:i]
	}
	if bytes.IndexByte(bc, 'N') != -1 {
		return true
	}
	if bytes.HasPrefix(bc, []byte("0_")) || bytes.HasSuffix(bc, []byte("_0")) || bytes.Contains(bc, []byte("_0_")) {
		return true
	}
	for _, p := range [4]string{"A00", "B00", "C00", "D00"} {
		if bytes.Contains(bc, []byte(p)) {
			return true
		}
	}
	return false
}

// BXTAG and VXTAG are the haplotagging-format barcode/validity tag
// prefixes, used by Haplotag2Corefq to locate a BX:Z: and an optional,
// explicit VX:i: field.
var BXTAG = []byte{'B', 'X', ':', 'Z', ':'}
var VXTAG = []byte{'V', 'X', ':', 'i', ':'}
