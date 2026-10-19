"""Random linked-read FASTQ for bx_stats_fq: usage: gen_fq.py PLATFORM SEED N OUT.fq

PLATFORM is haplotagging, stlfr or tellseq. Every read has a description (the Python original crashes on a
haplotagging read with none), and every barcode is well formed, so Python and djinn agree on the output.
Names carry an old-style /1 or /2 on some reads and a new-style "1:N:0:ATCG" description on others, so many
reads have both."""
import random, sys

platform, seed, n, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
rng = random.Random(seed)

def hap_bc(p_bad):
    return "".join(f"{l}{0 if rng.random() < p_bad else rng.randint(1, 96):02d}" for l in "ACBD")

def stlfr_bc(p_bad):
    return "_".join(str(0 if rng.random() < p_bad else rng.randint(1, 1500)) for _ in range(3))

def tellseq_bc(p_bad):
    return "".join("N" if rng.random() < p_bad else rng.choice("ACGT") for _ in range(18))

lines = []
for i in range(n):
    p_bad = rng.choice([0.0, 0.0, 0.02, 0.1, 0.3])   # a mix of clean and messy stretches
    base = f"LH00497:158:23LKF3LT4:{rng.randint(1, 8)}:{rng.randint(1101, 2400)}:{rng.randint(1, 50000)}:{rng.randint(1, 50000)}"
    has_bc = rng.random() < 0.9
    casava = rng.choice(["", " 1:N:0:ATGGAGA", " 2:Y:18:ACGTACGT"])
    if platform == "haplotagging":
        if has_bc:
            bc = hap_bc(p_bad)
            vx = 0 if "00" in bc else 1
            comment = rng.choice([f" BX:Z:{bc}", f" 1:N:0:ATGGAGA\tBX:Z:{bc}", f"\tVX:i:{vx}\tBX:Z:{bc}", f" BX:Z:{bc}\tVX:i:{vx}", f" 1:N:0:ATGGAGA\tVX:i:{vx}\tBX:Z:{bc}"])
        else:
            comment = rng.choice([" 1:N:0:ATGGAGA", " xx:Z:1", " VX:i:0"])
        header = f"@{base}{rng.choice(['', '/1', '/2'])}{comment}"
    elif platform == "stlfr":
        name = f"V300{rng.randint(1, 9999)}L1C001R{i:07d}"
        if has_bc:
            name += "#" + stlfr_bc(p_bad)
        old_style = rng.choice(["", "/1", "/2"])   # an old-style /1 or /2 on the name, and/or a new-style description
        header = f"@{name}{old_style}{casava}"
    else:
        name = base
        if has_bc:
            name += ":" + tellseq_bc(p_bad)
        header = f"@{name}{rng.choice(['', '/1', '/2'])}{casava}"
    seq = "".join(rng.choice("ACGT") for _ in range(30))
    lines += [header, seq, "+", "I" * 30]
open(out, "w").write("\n".join(lines) + "\n")
