package barcodes

// Generator is implemented by all barcode generator types, so one can be
// selected at runtime from a CLI flag while still writing into a
// caller-owned buffer instead of allocating per call.
type Generator interface {
	// NextInto writes the next barcode into dst (len(dst) >= MaxLen())
	// and returns bytes written and whether a barcode was produced.
	NextInto(dst []byte) (n int, ok bool)

	// InvalidInto writes the sentinel barcode into dst, returns bytes written.
	//	InvalidInto(dst []byte) int

	GetInvalid() []byte
	MaxLen() int
	Close()
}
