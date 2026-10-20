package fastq

import "testing"

func TestIsInvalidBarcode(t *testing.T) {
	tests := []struct {
		name string
		in   string
		want bool
	}{
		{"contains N", "ATCGN", true},
		{"haplotag A00", "A00C12B34D56", true},
		{"haplotag B00", "A12C12B00D56", true},
		{"haplotag D00 at end", "A12C12B34D00", true},
		{"stlfr leading zero", "0_5_5", true},
		{"stlfr middle zero", "5_0_5", true},
		{"stlfr trailing zero", "5_5_0", true},
		{"all valid haplotag", "A12C34B56D78", false},
		{"all valid stlfr", "5_5_5", false},
		{"pure ATCG", "ATCG", false},
		{"valid haplotag with hyphen suffix", "A12C34B56D78-1", false},
		{"otherwise-valid stlfr, but suffix contains N", "5_5_5-NN", false},
		{"invalid haplotag (A00) despite hyphen suffix", "A00C34B56D78-1", true},
		{"hyphen suffix never masks a real zero segment before it", "5_0_5-1", true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := isInvalidBarcode([]byte(tt.in))
			if got != tt.want {
				t.Errorf("isInvalidBarcode(%q) = %v, want %v", tt.in, got, tt.want)
			}
		})
	}
}
