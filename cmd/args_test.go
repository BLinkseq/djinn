package cmd

import (
	"os"
	"path/filepath"
	"testing"
)

// setStringFlag sets a string flag on cmd, returning a func that restores it
// to its previous value.
func setStringFlag(t *testing.T, cmdFlags interface{ Set(string, string) error }, name, value string) {
	t.Helper()
	if err := cmdFlags.Set(name, value); err != nil {
		t.Fatalf("setting flag %s=%s: %v", name, value, err)
	}
}

func makeReadableFile(t *testing.T, dir, name string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte("data"), 0o644); err != nil {
		t.Fatalf("setup file: %v", err)
	}
	return p
}

// ── sortXamCmd.Args ──────────────────────────────────────────────────────────

func TestSortXamArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")

	t.Run("no args", func(t *testing.T) {
		if err := sortXamCmd.Args(sortXamCmd, []string{}); err == nil {
			t.Errorf("expected error with no args")
		}
	})
	t.Run("too many args", func(t *testing.T) {
		if err := sortXamCmd.Args(sortXamCmd, []string{valid, valid, valid}); err == nil {
			t.Errorf("expected error with too many args")
		}
	})
	t.Run("nonexistent input file", func(t *testing.T) {
		if err := sortXamCmd.Args(sortXamCmd, []string{filepath.Join(dir, "missing.bam")}); err == nil {
			t.Errorf("expected error for missing input file")
		}
	})
	t.Run("valid input file", func(t *testing.T) {
		if err := sortXamCmd.Args(sortXamCmd, []string{valid}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
}

// ── convertXamCmd.Args ───────────────────────────────────────────────────────

func TestConvertXamArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")
	mapPath := filepath.Join(dir, "map.tsv")

	reset := func() {
		convertXamCmd.Flags().Set("map", "")
	}

	t.Run("no args", func(t *testing.T) {
		defer reset()
		if err := convertXamCmd.Args(convertXamCmd, []string{}); err == nil {
			t.Errorf("expected error with no args")
		}
	})
	t.Run("invalid target format", func(t *testing.T) {
		defer reset()
		setStringFlag(t, convertXamCmd.Flags(), "map", mapPath)
		if err := convertXamCmd.Args(convertXamCmd, []string{"bogus-format", valid}); err == nil {
			t.Errorf("expected error for unrecognized target format")
		}
	})
	t.Run("missing --map flag is required", func(t *testing.T) {
		defer reset()
		if err := convertXamCmd.Args(convertXamCmd, []string{"10x", valid}); err == nil {
			t.Errorf("expected error when --map is not provided")
		}
	})
	t.Run("valid args, format case-insensitive", func(t *testing.T) {
		defer reset()
		setStringFlag(t, convertXamCmd.Flags(), "map", mapPath)
		for _, format := range []string{"10x", "HAPLOTAGGING", "TellSeq", "stlfr"} {
			if err := convertXamCmd.Args(convertXamCmd, []string{format, valid}); err != nil {
				t.Errorf("format %q: unexpected error: %v", format, err)
			}
		}
	})
	t.Run("nonexistent input file", func(t *testing.T) {
		defer reset()
		setStringFlag(t, convertXamCmd.Flags(), "map", mapPath)
		if err := convertXamCmd.Args(convertXamCmd, []string{"10x", filepath.Join(dir, "missing.bam")}); err == nil {
			t.Errorf("expected error for missing input file")
		}
	})
}

// ── invalidXamCmd.Args / invalidFqCmd.Args ──────────────────────────────────

func TestInvalidXamArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")
	reset := func() { invalidXamCmd.Flags().Set("invalid", "") }

	t.Run("valid, no invalid flag", func(t *testing.T) {
		defer reset()
		if err := invalidXamCmd.Args(invalidXamCmd, []string{valid}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
	t.Run("invalid flag pointing at unwritable dir", func(t *testing.T) {
		defer reset()
		blocker := filepath.Join(dir, "blocker.txt")
		if err := os.WriteFile(blocker, []byte("x"), 0o644); err != nil {
			t.Fatalf("setup: %v", err)
		}
		setStringFlag(t, invalidXamCmd.Flags(), "invalid", filepath.Join(blocker, "sub", "prefix"))
		if err := invalidXamCmd.Args(invalidXamCmd, []string{valid}); err == nil {
			t.Errorf("expected error for unwritable --invalid destination")
		}
	})
}

func TestInvalidFqArgsRejectsSamePrefix(t *testing.T) {
	dir := t.TempDir()
	r1 := makeReadableFile(t, dir, "r1.fq")
	r2 := makeReadableFile(t, dir, "r2.fq")
	prefix := filepath.Join(dir, "out")
	reset := func() { invalidFqCmd.Flags().Set("invalid", "") }

	t.Run("same prefix for valid and invalid output is rejected", func(t *testing.T) {
		defer reset()
		setStringFlag(t, invalidFqCmd.Flags(), "invalid", prefix)
		if err := invalidFqCmd.Args(invalidFqCmd, []string{prefix, r1, r2}); err == nil {
			t.Errorf("expected error when --invalid prefix equals the valid-output prefix")
		}
	})
	t.Run("distinct prefixes accepted", func(t *testing.T) {
		defer reset()
		setStringFlag(t, invalidFqCmd.Flags(), "invalid", filepath.Join(dir, "bad"))
		if err := invalidFqCmd.Args(invalidFqCmd, []string{prefix, r1, r2}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
}

// ── unlinkedXamCmd.Args / unlinkedFqCmd.Args (singletons) ──────────────────

func TestUnlinkedXamArgsDefaultSingletonsIsOptional(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")
	reset := func() {
		unlinkedXamCmd.Flags().Set("singletons", "")
		unlinkedXamCmd.Flags().Set("bc-count", "")
	}
	defer reset()

	// Regression: --singletons is documented as optional (default ""), and
	// the Args func should not fail (or need write access at "."/cwd) when
	// it's left unset.
	if err := unlinkedXamCmd.Args(unlinkedXamCmd, []string{valid}); err != nil {
		t.Errorf("unexpected error with default (empty) --singletons: %v", err)
	}
}

func TestUnlinkedXamArgsWithSingletonsSet(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")
	reset := func() {
		unlinkedXamCmd.Flags().Set("singletons", "")
		unlinkedXamCmd.Flags().Set("bc-count", "")
	}
	defer reset()

	setStringFlag(t, unlinkedXamCmd.Flags(), "singletons", filepath.Join(dir, "singles.bam"))
	if err := unlinkedXamCmd.Args(unlinkedXamCmd, []string{valid}); err != nil {
		t.Errorf("unexpected error: %v", err)
	}
}

func TestUnlinkedFqArgsRejectsSamePrefix(t *testing.T) {
	dir := t.TempDir()
	r1 := makeReadableFile(t, dir, "r1.fq")
	r2 := makeReadableFile(t, dir, "r2.fq")
	prefix := filepath.Join(dir, "out")
	reset := func() {
		unlinkedFqCmd.Flags().Set("singletons", "")
		unlinkedFqCmd.Flags().Set("bc-count", "")
	}
	defer reset()

	setStringFlag(t, unlinkedFqCmd.Flags(), "singletons", prefix)
	if err := unlinkedFqCmd.Args(unlinkedFqCmd, []string{prefix, r1, r2}); err == nil {
		t.Errorf("expected error when --singletons prefix equals the valid-output prefix")
	}
}

// ── ncbiXam2FqCmd.Args ───────────────────────────────────────────────────────

func TestNcbiArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")

	if err := ncbiXam2FqCmd.Args(ncbiXam2FqCmd, []string{}); err == nil {
		t.Errorf("expected error with no args")
	}
	if err := ncbiXam2FqCmd.Args(ncbiXam2FqCmd, []string{"prefix"}); err == nil {
		t.Errorf("expected error with only one arg (needs prefix + file)")
	}
	if err := ncbiXam2FqCmd.Args(ncbiXam2FqCmd, []string{"prefix", filepath.Join(dir, "missing.bam")}); err == nil {
		t.Errorf("expected error for missing input file")
	}
	if err := ncbiXam2FqCmd.Args(ncbiXam2FqCmd, []string{"prefix", valid}); err != nil {
		t.Errorf("unexpected error: %v", err)
	}
}

// ── extractXamCmd.Args / extractFqCmd.Args ──────────────────────────────────

func TestExtractXamArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")

	if err := extractXamCmd.Args(extractXamCmd, []string{}); err == nil {
		t.Errorf("expected error with no args")
	}
	if err := extractXamCmd.Args(extractXamCmd, []string{valid, valid}); err == nil {
		t.Errorf("expected error with more than one arg (ExactArgs(1))")
	}
	if err := extractXamCmd.Args(extractXamCmd, []string{valid}); err != nil {
		t.Errorf("unexpected error: %v", err)
	}
}

func TestExtractFqArgs(t *testing.T) {
	dir := t.TempDir()
	r1 := makeReadableFile(t, dir, "r1.fq")
	r2 := makeReadableFile(t, dir, "r2.fq")

	if err := extractFqCmd.Args(extractFqCmd, []string{}); err == nil {
		t.Errorf("expected error with no args")
	}
	if err := extractFqCmd.Args(extractFqCmd, []string{r1, r2}); err != nil {
		t.Errorf("unexpected error: %v", err)
	}
}

// ── countXamCmd / countFqCmd ─────────────────────────────────────────────────

func TestCountXamArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")

	if err := countXamCmd.Args(countXamCmd, []string{}); err == nil {
		t.Errorf("expected error with no args")
	}
	if err := countXamCmd.Args(countXamCmd, []string{valid}); err != nil {
		t.Errorf("unexpected error: %v", err)
	}
}

// TestCountXamHasThreadsFlag is a regression test: countXamCmd's RunE reads
// the "threads" flag via cmd.Flags().GetInt("threads"), but until fixed, no
// "threads" flag was ever registered for countXamCmd in init(), so every
// invocation of `djinn sam count` failed immediately with
// "flag accessed but not defined: threads" before doing any work.
func TestCountXamHasThreadsFlag(t *testing.T) {
	if countXamCmd.Flags().Lookup("threads") == nil {
		t.Fatalf("countXamCmd is missing a registered --threads flag")
	}
	if _, err := countXamCmd.Flags().GetInt("threads"); err != nil {
		t.Errorf("GetInt(\"threads\") returned an error, RunE would fail: %v", err)
	}
}

// ── sampleXamCmd / sampleFqCmd downsample validation (in RunE) ─────────────

func TestSampleXamArgsRejectsNonPositiveDownsample(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")
	reset := func() {
		sampleXamCmd.Flags().Set("downsample", "0")
		sampleXamCmd.Flags().Set("seed", "-1")
	}
	defer reset()

	// The downsample <= 0 check lives in Args, not RunE -- cobra always
	// runs Args before RunE, but calling RunE directly (as this test
	// used to) skips Args entirely and ends up actually invoking
	// sample.SampleXam on the fake input, which doesn't itself validate
	// downsample and goes on to touch the filesystem for no reason this
	// test cares about.
	setStringFlag(t, sampleXamCmd.Flags(), "downsample", "0")
	if err := sampleXamCmd.Args(sampleXamCmd, []string{valid}); err == nil {
		t.Fatalf("expected error for downsample <= 0")
	}
}

func TestSampleFqArgsRejectsNonPositiveDownsample(t *testing.T) {
	dir := t.TempDir()
	// Args checks that r1/r2 exist before it checks downsample; use real
	// (if fake-content) files so the downsample check is what actually
	// gets exercised here, not a file-not-found error for an unrelated
	// reason -- see TestSampleXamArgsRejectsNonPositiveDownsample.
	r1 := makeReadableFile(t, dir, "r1.fq")
	r2 := makeReadableFile(t, dir, "r2.fq")
	reset := func() {
		sampleFqCmd.Flags().Set("downsample", "0")
	}
	defer reset()

	setStringFlag(t, sampleFqCmd.Flags(), "downsample", "-5")
	if err := sampleFqCmd.Args(sampleFqCmd, []string{"prefix", r1, r2}); err == nil {
		t.Fatalf("expected error for negative downsample value")
	}
}
