package cmd

import (
	"path/filepath"
	"testing"
)

func TestStatsXamArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")

	t.Run("no args", func(t *testing.T) {
		if err := statsXamCmd.Args(statsXamCmd, []string{}); err == nil {
			t.Errorf("expected error with no args")
		}
	})
	t.Run("too many args", func(t *testing.T) {
		if err := statsXamCmd.Args(statsXamCmd, []string{valid, valid}); err == nil {
			t.Errorf("expected error with too many args")
		}
	})
	t.Run("nonexistent input file", func(t *testing.T) {
		if err := statsXamCmd.Args(statsXamCmd, []string{filepath.Join(dir, "missing.bam")}); err == nil {
			t.Errorf("expected error for missing input file")
		}
	})
	t.Run("valid input file", func(t *testing.T) {
		if err := statsXamCmd.Args(statsXamCmd, []string{valid}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
	t.Run("stdin", func(t *testing.T) {
		if err := statsXamCmd.Args(statsXamCmd, []string{"-"}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
	t.Run("negative distance threshold", func(t *testing.T) {
		setStringFlag(t, statsXamCmd.Flags(), "distance-threshold", "-1")
		defer setStringFlag(t, statsXamCmd.Flags(), "distance-threshold", "0")
		if err := statsXamCmd.Args(statsXamCmd, []string{valid}); err == nil {
			t.Errorf("expected error for a negative distance threshold")
		}
	})
	t.Run("positive distance threshold", func(t *testing.T) {
		setStringFlag(t, statsXamCmd.Flags(), "distance-threshold", "50000")
		defer setStringFlag(t, statsXamCmd.Flags(), "distance-threshold", "0")
		if err := statsXamCmd.Args(statsXamCmd, []string{valid}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
}

func TestDeconvolveXamArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.bam")
	setDistance := func(t *testing.T, v string) {
		t.Helper()
		setStringFlag(t, deconvolveXamCmd.Flags(), "distance-threshold", v)
		t.Cleanup(func() { setStringFlag(t, deconvolveXamCmd.Flags(), "distance-threshold", "0") })
	}

	t.Run("no args", func(t *testing.T) {
		setDistance(t, "100")
		if err := deconvolveXamCmd.Args(deconvolveXamCmd, []string{}); err == nil {
			t.Errorf("expected error with no args")
		}
	})
	t.Run("too many args", func(t *testing.T) {
		setDistance(t, "100")
		if err := deconvolveXamCmd.Args(deconvolveXamCmd, []string{valid, valid}); err == nil {
			t.Errorf("expected error with too many args")
		}
	})
	t.Run("nonexistent input file", func(t *testing.T) {
		setDistance(t, "100")
		if err := deconvolveXamCmd.Args(deconvolveXamCmd, []string{filepath.Join(dir, "missing.bam")}); err == nil {
			t.Errorf("expected error for missing input file")
		}
	})
	t.Run("distance is required", func(t *testing.T) {
		// the flag defaults to 0, which would split nothing
		if err := deconvolveXamCmd.Args(deconvolveXamCmd, []string{valid}); err == nil {
			t.Errorf("expected error when no distance threshold was given")
		}
	})
	t.Run("negative distance", func(t *testing.T) {
		setDistance(t, "-5")
		if err := deconvolveXamCmd.Args(deconvolveXamCmd, []string{valid}); err == nil {
			t.Errorf("expected error for a negative distance")
		}
	})
	t.Run("valid", func(t *testing.T) {
		setDistance(t, "50000")
		if err := deconvolveXamCmd.Args(deconvolveXamCmd, []string{valid}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
	t.Run("stdin", func(t *testing.T) {
		setDistance(t, "50000")
		if err := deconvolveXamCmd.Args(deconvolveXamCmd, []string{"-"}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
}

func TestStatsFqArgs(t *testing.T) {
	dir := t.TempDir()
	valid := makeReadableFile(t, dir, "in.fq")

	t.Run("no args", func(t *testing.T) {
		if err := statsFqCmd.Args(statsFqCmd, []string{}); err == nil {
			t.Errorf("expected error with no args")
		}
	})
	t.Run("too many args", func(t *testing.T) {
		if err := statsFqCmd.Args(statsFqCmd, []string{valid, valid}); err == nil {
			t.Errorf("expected error with too many args")
		}
	})
	t.Run("nonexistent input file", func(t *testing.T) {
		if err := statsFqCmd.Args(statsFqCmd, []string{filepath.Join(dir, "missing.fq")}); err == nil {
			t.Errorf("expected error for missing input file")
		}
	})
	t.Run("valid input file, platform detected", func(t *testing.T) {
		if err := statsFqCmd.Args(statsFqCmd, []string{valid}); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
	for _, platform := range []string{"haplotagging", "stlfr", "tellseq", "StLFR", "TELLSEQ"} {
		t.Run("platform "+platform, func(t *testing.T) {
			setStringFlag(t, statsFqCmd.Flags(), "platform", platform)
			defer setStringFlag(t, statsFqCmd.Flags(), "platform", "")
			if err := statsFqCmd.Args(statsFqCmd, []string{valid}); err != nil {
				t.Errorf("unexpected error: %v", err)
			}
		})
	}
	t.Run("unknown platform", func(t *testing.T) {
		setStringFlag(t, statsFqCmd.Flags(), "platform", "10x")
		defer setStringFlag(t, statsFqCmd.Flags(), "platform", "")
		if err := statsFqCmd.Args(statsFqCmd, []string{valid}); err == nil {
			t.Errorf("expected error for an unknown platform")
		}
	})
}
