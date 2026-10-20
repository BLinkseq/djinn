---
title: Convert between FASTQ formats
sidebar:
  order: 4
  label: convert
---

import { Badge } from '@astrojs/starlight/components';



In the event you need your linked-read data converted into a different linked-read format, djinn has you covered.
We might disagree on the fragmented format landscape, but that doesn't mean you
shouldn't be able to use your data how and where you want to. This command converts a paired-end read set of FASTQ
files between the common linked-read FASTQ types.

## Usage
Accepts one unpaired FASTQ file or two paired FASTQ files (R1 and R2, can be gzipped), with barcodes specified in BX:Z SAM tag. Writes prefix.R1.fq.gz (and prefix.R2.fq.gz if paired). The first positional argument outFormat specifies the target data format (standard, haplotagging, tellseq, stlfr, 10x). Specifying 10X output will draw barcodes from the official 10X Genomics barcode whitelist.
```bash
djinn fastq convert PREFIX TARGET FQ1 FQ2
```

```bash example | tellseq → stlfr
djinn fastq convert data/orcs_stlfr stlfr data/orcs.R1.fq.gz data/orcs.R2.fq.gz
```

Auto-detects the input format as one of 10X, haplotagging, TELLseq, or stLFR,
and converts it to the format provided as the `TARGET` positional argument. If the
input data is in 10X format, where the barcode is the first [usually] 16 bases in the R1 sequence,
you will need to provide a `--barcodes` file that lets Djinn know what barcodes to look for.
In all cases, a file will be created with the barcode conversion map.

<details>
    <summary> example `--barcodes` file </summary>
```
ATGGAAGCCGTAGTTA
ACGGAAGCCGTAGTTC
ATGGAAGAAATAGTTA
ATGTTTGCCGTAGTTT
```
</details>


## Conversion targets

| `TARGET` | barcode format                                   | example      |
| :---------------------------- | :----------------------------------------------- | :------------------------------- |
| `10x`                         | the first N base pairs of R1 |                                  |
| `haplotagging`                | a `BX:Z:ACBD` SAM tag in the sequence header     | `@SEQID BX:Z:A01C93B56D11`       |
| `stlfr`                       | `#1_2_3` format appended to the sequence ID      | `@SEQID#1_2_3`                   |
| `tellseq`                     | `:ATCG` format appended to the sequence ID       | `@SEQID:GGCAAATATCGAGAAGTC`      |



## Options
| `-b` `--bc-map`  | Valid barcode list <Badge text="required: for 10X Genomics input" variant="tip" />   |
