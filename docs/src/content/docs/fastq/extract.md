---
title: Extract linked-read barcodes
sidebar:
  order: 4
  label: extract
---

You might want a list of all the unique barcodes in your data and `extract` will do just that.

### FASTQ
```bash usage
djinn fastq extract INPUT...
```

```bash example 
djinn fastq extract sample1.F.fq.gz sample1.R.fq.gz > sample1.bc
```
