package stats

import (
	"strconv"
	"strings"

	"github.com/biogo/hts/sam"
)

// ── record lengths ────────────────────────────────────────────────────────────
//
// These follow pysam's definitions (AlignedSegment.query_alignment_length,
// infer_query_length and reference_end) exactly, quirks included, because the
// molecule table is compared against output that was produced with pysam.

// alignedLength is the number of query bases between the soft clips
// (pysam's query_alignment_length). The query length comes from SEQ, or
// from the CIGAR when SEQ is absent.
func alignedLength(r *sam.Record) int {
	qlen := r.Seq.Length
	start := 0
	for _, op := range r.Cigar {
		if op.Type() != sam.CigarSoftClipped && op.Type() != sam.CigarHardClipped {
			break
		}
		if op.Type() == sam.CigarSoftClipped {
			start += op.Len()
		}
	}

	end := qlen
	if qlen == 0 {
		// no SEQ: pysam sums the CIGAR, counting a soft clip only if it comes before any aligned base
		for _, op := range r.Cigar {
			switch op.Type() {
			case sam.CigarMatch, sam.CigarInsertion, sam.CigarEqual, sam.CigarMismatch:
				end += op.Len()
			case sam.CigarSoftClipped:
				if end == 0 {
					end += op.Len()
				}
			}
		}
	} else {
		for i := len(r.Cigar) - 1; i >= 1; i-- { // pysam never looks at the first operation here
			op := r.Cigar[i]
			if op.Type() == sam.CigarHardClipped {
				continue
			}
			if op.Type() != sam.CigarSoftClipped {
				break
			}
			end -= op.Len()
		}
	}
	return end - start
}

// inferredQueryLength is the read length implied by the CIGAR, without
// hard-clipped bases (pysam's infer_query_length).
func inferredQueryLength(r *sam.Record) int {
	n := 0
	for _, op := range r.Cigar {
		switch op.Type() {
		case sam.CigarMatch, sam.CigarInsertion, sam.CigarSoftClipped, sam.CigarEqual, sam.CigarMismatch:
			n += op.Len()
		}
	}
	return n
}

// alignedEnd is one past the last reference base covered (pysam's
// reference_end, i.e. htslib's bam_endpos): pos+1 when the CIGAR consumes no reference.
func alignedEnd(r *sam.Record) int {
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

// insertSize is the insert length a read contributes to its molecule.
func insertSize(r *sam.Record) int {
	if r.Flags&sam.Paired != 0 {
		if r.Flags&sam.Supplementary != 0 {
			return alignedLength(r)
		}
		return max(0, r.TempLen)
	}
	return max(abs(r.TempLen), inferredQueryLength(r))
}

func abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

// ── number formatting ─────────────────────────────────────────────────────────

// pyFloat formats f like Python's repr(float): shortest round-trip digits,
// ".0" on whole numbers, and an exponent below 1e-4 ("3e-05").
func pyFloat(f float64) string {
	s := strconv.FormatFloat(f, 'g', -1, 64)
	if !strings.ContainsAny(s, ".e") {
		s += ".0"
	}
	return s
}

// coverage is max(0, round(min(num/denom, 1), 5)), rounded the way Python's
// round() does (correctly, from the exact binary value).
func coverage(num, denom int) float64 {
	v := min(float64(num)/float64(denom), 1.0)
	v, _ = strconv.ParseFloat(strconv.FormatFloat(v, 'f', 5, 64), 64)
	if !(v > 0) {
		return 0
	}
	return v
}
