#!/bin/bash
# Regenerates the test inputs and their expected outputs from the Python originals
# (harpy/utils/bx_stats_sam.py and bx_stats_fq.py; needs pysam and click):
#
#   HARPY_UTILS=/path/to/harpy/harpy/utils ./regen.sh
#
# The expected files are the specification: they come from the Python code, not from this package.
set -euo pipefail
: "${HARPY_UTILS:?set HARPY_UTILS to the harpy/harpy/utils directory}"
cd "$(dirname "$0")"
export LC_ALL=C
python3 gen_sam.py 5 110 small.sam
for d in 0 1 300 2000 40000; do
    python3 - "$HARPY_UTILS" "$d" <<'PY' | sort > "small.d$d.tsv"
import sys
sys.path.insert(0, sys.argv[1])
from bx_stats_sam import bx_stats_sam
bx_stats_sam(args=["-d", sys.argv[2], "small.sam"], prog_name="bx_stats_sam")
PY
done

# FASTQ: one small file per platform, and what bx_stats_fq.py reports for it
mkdir -p fq
for platform in haplotagging stlfr tellseq; do
    python3 gen_fq.py "$platform" 5 250 "fq/$platform.fq"
    python3 oracle_fq.py "$HARPY_UTILS" "$platform" "fq/$platform.fq" > "fq/$platform.tsv"
done
