---
title: Count linked-read barcodes
sidebar:
  order: 4
  label: count
---

You might want a table of the counts (frequency) of all the unique barcodes in your data and that's what [!badge corners="pill" text="count"] is for. Records
with missing barcodes are ignored.

Counting barcodes in SAM/BAM files avoids double-counting by counting forward reads if they are unpaired or unpaired, and reverse reads if they are unpaired or the barcode has not been observed yet.
```bash usage
djinn sam count INPUT
```

```bash example
djinn sam count sample1.bam > sample1.bxcount
```
