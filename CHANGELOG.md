## New
- `sam mv-bx` to move `BX:Z` barcodes to the end of the record 
- `fq mv-bx` to move `BX:Z` barcodes to the end of the read header
 
## Changed
- commands that write new files where barcodes are modified (e.g., `convert`) put the `BX:Z` tag at the end (if applicable)
