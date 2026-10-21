Rewritten with Go for speed and resource efficiency.

## Changes
Just about everything. New CLI. It's _similar_, but not the same. **Please read the docs!!!**

## Breaking
Almost everything

## New
- linked-read stats provided by Harpy moved here
  - `djinn sam stats` and `djinn fq stats`
- comprehensive 10X support
  - 10X barcode whitelist embedded into the program, but one can be provided if it differs
- no external dependencies -- sorting is no longer a round-trip to `samtools`
- now also provides various Go packages for processessing FASTQ and SAM linked-read data in other applications
