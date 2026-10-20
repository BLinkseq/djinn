package sort

const chunkSize = 500_000 // records per temp file, tune to RAM
const ioBatchSize = 2000  // records per channel send on the read/write hot paths

// missingBXKey sorts after every real BX value (samtools convention: unset
// tags trail). 0xff can't appear in a real BX:Z value (SAM aux text is
// printable ASCII), so this can't collide with a real key.
var missingBXKey = []byte{0xff, 0xff, 0xff, 0xff}
