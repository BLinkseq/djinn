#!/bin/bash
# Regenerates small.sam and the expected outputs from the Python original
# (harpy/utils/bx_stats_sam.py; needs pysam and click):
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
