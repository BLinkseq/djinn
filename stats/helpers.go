package stats

import (
	"strconv"

	"github.com/biogo/hts/sam"
)

// ── measuring a read ──────────────────────────────────────────────────────────
//
// Alignment files describe how each read lines up with the reference genome using a short
// "recipe" called a CIGAR, made of operations like 50M (50 bases line up), 5S (5 bases at the
// end that were left out, "soft clipped"), 2I (2 extra bases in the read) or 3D (3 bases missing
// from the read). The functions below read that recipe to measure a read.
//
// They copy the exact rules of pysam (the library the original script used), quirks
// included, because the numbers must match the original program's output digit for digit.

// alignedLength is the number of bases of the read that actually line up with the reference,
// i.e. the read without its soft-clipped ends (pysam's query_alignment_length).
func alignedLength(r *sam.Record) int {
	readLen := r.Seq.Length // how many bases the record stores (0 if the sequence was left out)

	// Count the soft-clipped bases at the start of the read. Hard clips (H) are bases that were
	// cut out of the file entirely, so they take up no room in the read and are stepped over.
	start := 0
	for _, op := range r.Cigar {
		if op.Type() != sam.CigarSoftClipped && op.Type() != sam.CigarHardClipped {
			break // reached the first operation that is neither kind of clip
		}
		if op.Type() == sam.CigarSoftClipped {
			start += op.Len()
		}
	}

	end := readLen // where the aligned part of the read stops
	if readLen == 0 {
		// The sequence was left out of the file, so work the read's length out from the recipe.
		// pysam counts a soft clip here only if it comes before any aligned base (a leading clip),
		// so clips at the far end are not counted.
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
		// The sequence is present, so its length is known; trim the soft clips off the far end.
		// Walk the recipe backwards, skipping hard clips, until something that is not a clip.
		// (pysam never looks at the very first operation in this loop, hence i >= 1.)
		for i := len(r.Cigar) - 1; i >= 1; i-- {
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

// inferredQueryLength is the length of the read worked out from the recipe, counting every
// base the read contains but not hard-clipped ones, which are not part of it
// (pysam's infer_query_length).
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

// insertSize is the length of the original DNA fragment the read came from, as best we can tell.
// For paired reads the file states it (the "template length"), which is negative on the second
// mate, so we clamp at zero. For an unpaired read we take the larger of the stated length and the
// read's own length.
func insertSize(r *sam.Record) int {
	if r.Flags&sam.Paired != 0 {
		if r.Flags&sam.Supplementary != 0 {
			return alignedLength(r) // an extra piece of a read: all we know is how much of it lined up
		}
		return max(0, r.TempLen)
	}
	return max(abs(r.TempLen), inferredQueryLength(r))
}

// abs is the absolute value of an int (Go's math.Abs works only on floating-point numbers).
func abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

// appendCoverage adds num/denom to dst as a fraction with 5 decimal places, capped at 1.
func appendCoverage(dst []byte, num, denom int) []byte {
	return strconv.AppendFloat(dst, min(float64(num)/float64(denom), 1.0), 'f', 5, 64)
}
