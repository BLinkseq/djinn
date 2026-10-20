Rewritten with Go for speed and resource efficiency.

### SAM Benchmarks
There is a noticeable and somewhat expected performance
drop in `sort` due to being implemented natively in Go
and slower than the C implementation that was
in samtools. This was added to remove the samtools dep, but
you're welcome to continue using `samtools sort -t BX`.
Invalid filtering is a bit slower as well, and it's likely
due to the native Go compression algorithm in reading/writing
BAM. Everything else is noticebly improved.

| function | Python | Go | net change |
|:---------|-------:|-----:|:--------:|
| count | 6.85 | 1.50 | +356.6% |
| extract | 4.33 | 1.44 | +200.6% |
| filter-invalid | 2.38 | 3.34 | -28.6% |
| filter-singletons | incomplete | 7.67 | +1000%+ |
| ncbi | 2.55 | 2.37 | +7.6% |
| sample | 6.60 | 5.1 | +29.4% |
| sort | 3.63 | 9.27 | -63.7% |
| standardize | 12.97 | 3.80 | +241.3% |
