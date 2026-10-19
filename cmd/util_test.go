package cmd

import (
	"os"
	"path/filepath"
	"testing"
)

func TestFileExists(t *testing.T) {
	dir := t.TempDir()
	existing := filepath.Join(dir, "present.txt")
	if err := os.WriteFile(existing, []byte("x"), 0o644); err != nil {
		t.Fatalf("setup: %v", err)
	}
	missing := filepath.Join(dir, "missing.txt")

	cases := []struct {
		name    string
		path    string
		wantErr bool
	}{
		{"existing file", existing, false},
		{"missing file", missing, true},
		{"a directory", dir, true},
		// A path os.Stat cannot even parse (NUL byte) previously caused a
		// nil-pointer dereference: os.Stat returns a non-nil, non-NotExist
		// error with info == nil, and the old code unconditionally called
		// info.IsDir() regardless of that error.
		{"unstat-able path", "bad\x00path", true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			err := fileExists(c.path)
			if (err != nil) != c.wantErr {
				t.Errorf("fileExists(%q) error = %v, wantErr %v", c.path, err, c.wantErr)
			}
		})
	}
}

func TestFileReadable(t *testing.T) {
	dir := t.TempDir()
	existing := filepath.Join(dir, "present.txt")
	if err := os.WriteFile(existing, []byte("x"), 0o644); err != nil {
		t.Fatalf("setup: %v", err)
	}
	missing := filepath.Join(dir, "missing.txt")

	if err := fileReadable(existing); err != nil {
		t.Errorf("fileReadable(existing) = %v, want nil", err)
	}
	if err := fileReadable(missing); err == nil {
		t.Errorf("fileReadable(missing) = nil, want error")
	}
}

func TestFilecheck(t *testing.T) {
	dir := t.TempDir()
	existing := filepath.Join(dir, "present.txt")
	if err := os.WriteFile(existing, []byte("x"), 0o644); err != nil {
		t.Fatalf("setup: %v", err)
	}
	missing := filepath.Join(dir, "missing.txt")

	if err := filecheck(existing); err != nil {
		t.Errorf("filecheck(existing) = %v, want nil", err)
	}
	if err := filecheck(missing); err == nil {
		t.Errorf("filecheck(missing) = nil, want error")
	}
	if err := filecheck(dir); err == nil {
		t.Errorf("filecheck(dir) = nil, want error (directory)")
	}
}

func TestCheckIfExecInPath(t *testing.T) {
	if !checkIfExecInPath("ls") {
		t.Errorf("expected 'ls' to be found in PATH")
	}
	if checkIfExecInPath("this-executable-should-not-exist-anywhere-xyz") {
		t.Errorf("expected bogus executable to not be found in PATH")
	}
}

func TestEnsureWritableDir(t *testing.T) {
	dir := t.TempDir()

	t.Run("writable existing dir", func(t *testing.T) {
		target := filepath.Join(dir, "sub", "file.txt")
		if err := ensureWritableDir(target); err != nil {
			t.Errorf("ensureWritableDir(%q) = %v, want nil", target, err)
		}
		if info, err := os.Stat(filepath.Join(dir, "sub")); err != nil || !info.IsDir() {
			t.Errorf("expected parent directory to be created")
		}
	})

	t.Run("directory not creatable under a file", func(t *testing.T) {
		// dir/blocker is a regular file; using it as a parent directory
		// component must fail.
		blocker := filepath.Join(dir, "blocker")
		if err := os.WriteFile(blocker, []byte("x"), 0o644); err != nil {
			t.Fatalf("setup: %v", err)
		}
		target := filepath.Join(blocker, "sub", "file.txt")
		if err := ensureWritableDir(target); err == nil {
			t.Errorf("ensureWritableDir(%q) = nil, want error", target)
		}
	})
}
