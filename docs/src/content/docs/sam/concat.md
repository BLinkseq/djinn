---
title: Concatenate linked-read SAM files
sidebar:
  order: 4
  label: concat
---
  

If you want to combine linked-read data, you run into an immediate issue: barcode clashing. Different 
samples will likely share barcodes by chance. A single sample that was made into multiple linked-read
libraries will have barcode clashing across libraries for the same reason. Therefore, concatenating needs
to happen with barcode recoding to make sure barcodes are completely unique between all things being merged.

```bash usage
djinn sam concat INPUT...  > output.bam
```

Concatenate records from linked-read SAM/BAM files while making sure molecule identification tags (`MI` or `BX`)
remain unique for every sample. This is a means of accomplishing the same as 'samtools cat', except all MI/BX tags are updated
so individuals don't have overlapping tags (which would mess up all the linked-read info). The default ignores existing `MI` tags
and writes new ones that correspond to unique `BX` tags. Using `--mi` is the opposite, where it ignores existing `BX` tags and writes
new ones that correspond with the `MI` tags in the barcode style of your choice.

```bash example (single sample, multiple libraries)
djinn sam concat sample1.lib1.bam sample1.lib2.bam sample1.lib3.bam > sample1.bam
```

```bash example (multi-sample)
djinn sam concat sample*.bam > pop_A.bam
```
