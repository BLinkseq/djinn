// Package deconvolve separates reads that share a barcode but came from different molecules.
//
// # Plain-language picture of what happens here
//
// A linked-read barcode is meant to say "these reads came from the same long DNA molecule".
// But a library holds many molecules, so by chance two different molecules can carry the same
// barcode. They can be told apart by where they sit on the chromosome: reads from one molecule
// huddle together, and a long empty stretch before the next read with that barcode suggests it
// came from a different molecule.
//
// This package walks through an alignment file from the start of each chromosome to its end. For
// every barcode it remembers where its latest read ended. When the next read with that barcode starts
// more than a chosen distance further along, it is taken to be a new molecule, and from then on the
// barcode is written with a number on the end: ACGT, then ACGT-1, then ACGT-2, and so on. Everything
// else about every record is left exactly as it was.
package deconvolve

import (
	"errors"
	"os"
	"strconv"

	"github.com/blinkseq/djinn/xam"

	"github.com/biogo/hts/sam"
)

// ErrNotSorted is returned when the input is not coordinate sorted.
// "Coordinate sorted" means the reads are listed in order of where they sit on each
// chromosome. The whole method relies on that order, so we refuse to guess without it.
var ErrNotSorted = errors.New("input alignment file does not appear to be coordinate-sorted, " +
	"which will result in incorrect output. Please sort the alignments by coordinates (e.g. samtools sort) and try again")

// DeconvolveXam reads a coordinate-sorted SAM/BAM file (or stdin for "-") and writes every record to
// stdout, with the barcode (BX tag) suffixed by "-1", "-2", ... wherever a barcode's reads are more
// than distance bases apart. The first molecule of a barcode keeps the barcode as it is. Barcode
// numbering starts over on each chromosome. Records without a valid barcode, and unmapped records,
// are written unchanged.
//
// A summary goes to stderr when the file is done: the number of molecules before and after, and
// how many barcodes were split how many times. A molecule here is one barcode on one chromosome.
//
// asSam writes SAM text instead of BAM. distance must be greater than 0.
func DeconvolveXam(infile string, distance int, asSam bool, threads int) error {
	if distance <= 0 {
		return errors.New("the distance threshold must be greater than 0")
	}
	readThread, writeThread := xam.IoThreads(threads)

	// A background worker reads the file and hands us records through a channel (a conveyor belt).
	recChan, br := xam.NewXamReaderChan(infile, xam.ChanCap, xam.IoBuf, readThread)

	// Note in the file's header that djinn touched it, as the other djinn commands do.
	hdr := br.Header()
	if err := hdr.AddProgram(xam.NewPG(hdr, "djinn sam deconvolve -d "+strconv.Itoa(distance)+" "+infile)); err != nil {
		return err
	}

	// Another background worker takes finished records from a second channel and writes them out.
	writeChan, writeDone := xam.NewXamWriterChan("-", hdr, xam.ChanCap, xam.IoBuf, writeThread, asSam)

	// The splitter says which molecule each record is in. At the end of each chromosome we count its
	// barcodes (and how often each was split) for the summary.
	var rep report
	split := NewSplitter(distance)
	split.OnContigEnd(func() {
		split.Barcodes(func(_ string, splits int) { rep.add(splits) })
	})

	var runErr error
	for rec := range recChan {
		place, err := split.Place(rec)
		if err != nil {
			runErr = err
			break
		}
		if place.Barcode != "" { // records without a usable barcode have nothing to suffix
			setSuffixed(rec, place.Barcode, place.Molecule)
		}
		writeChan <- rec // every record is written, whether or not its barcode changed
	}
	close(writeChan)
	<-writeDone // wait until everything has really been written out

	if runErr != nil {
		// We stopped early, so the reader may still be holding unread records. Let it finish up in the
		// background instead of leaving it stuck waiting for us to take the next one.
		go func() {
			for range recChan {
			}
		}()
		return runErr
	}
	split.Finish() // the last chromosome is still open
	return rep.write(os.Stderr)
}

// setSuffixed writes the barcode with "-n" on the end into the record's BX tag. Number 0 is the barcode's
// first molecule, which keeps the barcode exactly as it was (the record is left untouched).
func setSuffixed(rec *sam.Record, bx string, n int) {
	if n == 0 {
		return
	}
	suffixed := bx + "-" + strconv.Itoa(n)
	xam.SetBxString(rec, &suffixed) // replaces the BX tag, or adds one if the barcode came from the read name
}
