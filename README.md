![djinn logo](https://raw.githubusercontent.com/pdimens/djinn/refs/heads/docs/static/djinn.png)

[![GitHub Release](https://img.shields.io/github/v/release/pdimens/djinn?style=for-the-badge&logo=anaconda&logoColor=ffffff)](https://github.com/pdimens/djinn/releases)
[![documentation badge](https://img.shields.io/badge/read%20the-docs-fbab3a?style=for-the-badge&logo=quicklook&logoColor=ffffff)](https://pdimens.github.io/djinn)
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pdimens/djinn/tests.yaml?style=for-the-badge&logo=cachet&logoColor=ffffff)](https://www.youtube.com/watch?v=F1qdBPlK9M4)

Modify and convert linked-read FASTQ and BAM files

## Extract barcodes
Get all the barcodes from the input file(s)

## Sort by barcode
Sort records by barcode rather than name or position

## Molecule statistics
`djinn sam stats` summarizes each molecule of a coordinate-sorted SAM/BAM file: number of reads, start, end,
inferred length, aligned bp, inferred insert length and coverage. With `-d DISTANCE`, reads that share a barcode
are split into separate molecules (`BARCODE`, `BARCODE-1`, ...) when they are more than `DISTANCE` bases apart,
using the same splitting as `djinn sam deconvolve`. A barcode is valid when its `VX` tag says so, or, when there is
no `VX` tag, when the barcode itself is valid. `djinn fq stats` does the same kind of summary for a FASTQ file:
reads, reads with a barcode, valid and invalid barcodes, and how many barcodes were bad at each position
(haplotagging segment `A00`, stLFR part `0`, TELLseq base `N`).

## Deconvolution
Different molecules can share a barcode. `djinn sam deconvolve -d DISTANCE` tells them apart by where their reads
sit: when a barcode's reads are more than `DISTANCE` bases apart, the later ones get the barcode with a hyphenated
number on the end (`BARCODE`, `BARCODE-1`, `BARCODE-2`, ...). Everything else in every record is left alone, and a
summary (molecules before and after, and how many barcodes were split how many times) goes to stderr. Chain it into
`stats` to get per-molecule numbers: `djinn sam deconvolve -d 50000 in.bam | djinn sam stats -` (or just
`djinn sam stats -d 50000 in.bam`, which splits internally).

## Filtering
- separate linked reads from singletons
- omit invalid barcodes

## Convert files
Djinn converts between linked-read data formats. You can convert between formats in terms of FASTQ type or barcode style. It supports:
- 10X
- haplotagging
- stLFR
- TELLseq
- [standard](https://blinkseq.github.io/lastq/)

## NCBI submission
NCBI strips out sequence headers from FASTQ submissions, so it would be best to convert your linked-read
FASTQ data into an unaligned BAM file, with the linked-read barcode stored in the `BX` or `BC` tag.
Djinn provides a convenience function to convert to (or from) this format, although it's really just
a basic `samtools` command.

## Hi-C spoofing
This is **extremely experimental** and it converts a paired-end linked-read fastq file pair into one that conforms to Hi-C expectations by mix-matching the R1s and R2s of reads that share a barcode.


---
## Go Conversion Roadmap

| fastq | sam | function |
|:-----:|:-----:|:-----|
| [x] | [x] | convert           |     
| [x] | [x] | count             |     
| [x] | [x] | extract           |              
| [x] | [x] | filter-invalid    |            
| [x] | [x] | filter-singletons |   
| [x] | [x] | ncbi              | 
| [x] | [x] | sample            |               
| [x] | [x] | stats             |
| [ ] | [x] | sort              |             
| [x] | [x] | standardize       | 
| [ ] | [ ] | concat            |
| [ ] |     | spoof-hic         |               
|     | [x] | deconvolute       |
