package fastq

import (
	"fmt"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// parseHeader runs a parser on a single FASTQ record with the given header (everything after the "@").
func parseHeader(t *testing.T, parse func(*fastx.Record) (CoreFq, bool), header string) (CoreFq, bool) {
	t.Helper()
	seq.ValidateSeq = false
	path := writeFastqFile(t, "@"+header, "ACGT", "+", "IIII")
	r, err := fastx.NewReader(seq.DNA, path, "")
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	rec, err := r.Read()
	if err != nil {
		t.Fatal(err)
	}
	return parse(rec)
}

// Every combination of the three kinds of suffix a read can carry must be taken apart correctly:
//   - an old-style paired-end marker on the name ("/1", "/2"),
//   - a new-style marker in the description ("1:N:0:ATCG"),
//   - a hyphen suffix on the barcode itself ("-1", "-GEM1"), which belongs to the barcode.
//
// Whatever mix is present, the name comes out clean, the barcode keeps its own hyphen suffix and nothing
// else, and the marker is reported (the new-style one if both are there).
func TestParsers_SuffixCombinations(t *testing.T) {
	oldStyle := []struct{ suffix, marker string }{{"", ""}, {"/1", "1"}, {"/2", "2"}}
	newStyle := []string{"", "1:N:0:ATCG", "2:Y:18:ACGT"}
	hyphen := []string{"", "-1", "-GEM1"}

	for _, o := range oldStyle {
		for _, n := range newStyle {
			for _, h := range hyphen {
				wantCasava := o.marker
				if n != "" {
					wantCasava = n // the new-style marker wins when both are present
				}
				desc := ""
				if n != "" {
					desc = " " + n
				}

				for _, c := range []struct {
					name   string
					parse  func(*fastx.Record) (CoreFq, bool)
					header string
					bc     string
				}{
					{"stlfr", Stlfr2Corefq, "read#12_34_56" + h + o.suffix + desc, "12_34_56" + h},
					{"tellseq", Tellseq2Corefq, "read:ACGTACGT" + h + o.suffix + desc, "ACGTACGT" + h},
					// haplotagging: the tag is in the description, either before or after the marker
					{"haplotagging", Haplotag2Corefq, "read" + o.suffix + desc + "\tBX:Z:A01C02B03D04" + h + "\tVX:i:1", "A01C02B03D04" + h},
					{"haplotagging, tag first", Haplotag2Corefq, "read" + o.suffix + "\tBX:Z:A01C02B03D04" + h + desc, "A01C02B03D04" + h},
				} {
					label := fmt.Sprintf("%s %q", c.name, c.header)
					core, ok := parseHeader(t, c.parse, c.header)
					if !ok {
						t.Errorf("%s: no barcode found", label)
						continue
					}
					if string(core.ID) != "read" {
						t.Errorf("%s: id = %q, want %q", label, core.ID, "read")
					}
					if string(core.Barcode) != c.bc {
						t.Errorf("%s: barcode = %q, want %q", label, core.Barcode, c.bc)
					}
					if string(core.CASAVA) != wantCasava {
						t.Errorf("%s: casava = %q, want %q", label, core.CASAVA, wantCasava)
					}
					if len(core.Comments) != 0 {
						t.Errorf("%s: leftover description %q", label, core.Comments)
					}
				}
			}
		}
	}
}

// Whatever suffixes are present, technology detection must see the barcode.
func TestCoreFqParser_DetectsWithEverySuffixCombination(t *testing.T) {
	for _, suffix := range []string{"", "/1", "/2"} {
		for _, desc := range []string{"", " 1:N:0:ATCG"} {
			for _, hy := range []string{"", "-1"} {
				for want, header := range map[string]string{
					"stlfr":        "read#12_34_56" + hy + suffix + desc,
					"tellseq":      "read:ACGTACGT" + hy + suffix + desc,
					"haplotagging": "read" + suffix + desc + "\tBX:Z:A01C02B03D04" + hy,
				} {
					_, got, err := CoreFqParser(writeFastqFile(t, "@"+header, "ACGT", "+", "IIII"))
					if err != nil || got != want {
						t.Errorf("%q: detected %q (err %v), want %q", header, got, err, want)
					}
				}
			}
		}
	}
}

func TestSplitCASAVA(t *testing.T) {
	for _, c := range []struct {
		name, id, desc          string
		wantID, wantDesc, wantC string
	}{
		{"neither", "read", "BX:Z:AA", "read", "BX:Z:AA", ""},
		{"old style only", "read/1", "BX:Z:AA", "read", "BX:Z:AA", "1"},
		{"new style only", "read", "1:N:0:ATCG BX:Z:AA", "read", "BX:Z:AA", "1:N:0:ATCG"},
		{"both: both are removed, the new style is reported", "read/2", "1:N:0:ATCG\tBX:Z:AA", "read", "BX:Z:AA", "1:N:0:ATCG"},
		{"a slash that is not a marker stays", "read/3", "", "read/3", "", ""},
	} {
		id, desc, casava := splitCASAVA([]byte(c.id), []byte(c.desc))
		if string(id) != c.wantID || string(desc) != c.wantDesc || string(casava) != c.wantC {
			t.Errorf("%s: got (%q, %q, %q), want (%q, %q, %q)", c.name, id, desc, casava, c.wantID, c.wantDesc, c.wantC)
		}
	}
}
