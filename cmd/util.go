package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"

	"github.com/blinkseq/djinn/xam"
)

// check if the file exists and the directory is writable, return an
// error if not, otherwise nil
func filecheck(filename string) error {
	//var err error
	if err := fileExists(filename); err != nil {
		return err
	}
	if err := fileReadable(filename); err != nil {
		return err
	}
	return nil
}

// check if file exists. Returns a true if it does, otherwise false.
func fileExists(filename string) error {
	info, err := os.Stat(filename)
	if err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("\033[33;1m%s\033[0m does not exist", filename)
		}
		return err
	}
	if info.IsDir() {
		return fmt.Errorf("\033[33;1m%s\033[0m is a directory", filename)
	}
	return nil
}

// check if the input file has read permissions
func fileReadable(filename string) error {
	file, err := os.Open(filename)
	if err != nil {
		return fmt.Errorf("\033[33;1m%s\033[0m does not have read persmissions.\n", filename)
	}
	file.Close()
	return nil
}

// Check if the input path has write permissions
func ensureWritableDir(path string) error {
	dir := filepath.Dir(path)
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	f, err := os.CreateTemp(dir, ".writetest")
	if err != nil {
		return err
	}
	f.Close()
	return os.Remove(f.Name())
}

// Parse the thread count provided to `threads` at the CLI and clamp it between 1 and the system max
func safethreads(threads int) int {
	maxCores := runtime.NumCPU()
	// clamp between 1 and max system threads
	if threads > maxCores {
		fmt.Fprintf(
			os.Stderr,
			"Threads requested (%v) exceed total available (%v), scaling down to maximum available threads.\n",
			threads, maxCores,
		)
	}
	return min(maxCores, max(threads, 1))
}

// Check if the SAM/BAM file is coordinate sorted from its header SO value. Must match SO:coordinate,
// otherwise returns an error.
func isSamSorted(infile string) error {
	_, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, 1)
	// Note in the file's header that djinn touched it, as the other djinn commands do.
	if so := br.Header().SortOrder.String(); so != "coordinate" {
		return fmt.Errorf("Input must be coordinate sorted, but file header reports SO:%s", so)
	}
	return nil
}
