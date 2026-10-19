"""Random coordinate-sorted SAM exercising every branch of the molecule statistics: usage: gen_sam.py SEED N OUT.sam [OUT.bam]

A VX tag is only ever left out of reads that have no barcode, because the Python original treats a missing VX as
invalid while djinn infers validity from the barcode itself; with VX present on every barcoded read both agree."""
import random, sys, pysam

seed, n, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rng = random.Random(seed)
contigs = [("c1", 3_000_000), ("chr2", 2_000_000), ("3L", 1_000_000)]
bcs = [f"A{rng.randrange(100):02d}C{rng.randrange(100):02d}B{rng.randrange(100):02d}D{rng.randrange(100):02d}" for _ in range(60)]

def cigar():
    L, a, b, h, x, y, i, d, k = (rng.randint(1, 120), rng.randint(1, 30), rng.randint(1, 30), rng.randint(1, 20),
                                 rng.randint(1, 60), rng.randint(1, 60), rng.randint(1, 8), rng.randint(1, 15), rng.randint(1, 5))
    nn, p = rng.randint(50, 4000), rng.randint(1, 5)
    return rng.choice([
        (f"{L}M", L), (f"{a}S{L}M", a + L), (f"{L}M{b}S", L + b), (f"{a}S{L}M{b}S", a + L + b),
        (f"{x}M{i}I{y}M", x + i + y), (f"{x}M{d}D{y}M", x + y), (f"{x}M{nn}N{y}M", x + y),
        (f"{h}H{a}S{L}M", a + L), (f"{L}M{b}S{h}H", L + b), (f"{h}H{L}M", L), (f"{L}={k}X{x}=", L + k + x),
        (f"{x}M{p}P{y}M", x + y), (f"{x}M{d}D{i}I{y}M", x + i + y), (f"{a}S{x}M{d}D{y}M{b}S", a + x + y + b),
        (f"{a}I{L}M", a + L), (f"{L}M{d}D", L),
    ])

lines = ["@HD\tVN:1.6\tSO:coordinate"] + [f"@SQ\tSN:{c}\tLN:{l}" for c, l in contigs]
for cname, clen in contigs:
    pos = 1
    for r in range(n):
        pos += rng.choices([0, rng.randint(1, 300), rng.randint(300, 5000), rng.randint(5000, 60000)], [35, 50, 10, 5])[0]
        if pos > clen - 2000:
            break
        cg, ql = cigar()
        cat = rng.random()
        flag, cig = 0, cg
        if cat < 0.03:   flag, cig = 4, "*"                       # unmapped
        elif cat < 0.05: flag = 0x400                              # duplicate
        elif cat < 0.07: flag = 0x100                              # secondary
        elif cat < 0.09: cig = "*"                                 # mapped, no CIGAR
        elif cat < 0.17: flag = 0x800 | rng.choice([0, 1])         # supplementary (maybe paired)
        else:
            if rng.random() < 0.6:
                flag = 1 | 2 | rng.choice([0x40, 0x80]) | rng.choice([0, 16, 32])
            else:
                flag = rng.choice([0, 16])
        if flag & 1 and flag & 0x800:
            flag |= rng.choice([0x40, 0x80])
        if cig == "*" or rng.random() < 0.08:
            seq = qual = "*"
        else:
            seq = "".join(rng.choice("ACGT") for _ in range(ql)); qual = "I" * ql
        tlen = rng.choice([0, rng.randint(1, 800), -rng.randint(1, 800), rng.randint(0, 300)])
        tags = []
        e = rng.random()
        if e < 0.05: pass                                           # no BX at all
        else:
            bx = rng.choice(bcs)
            if rng.random() < 0.02: bx = rng.choice(bcs) + "_x"
            tags.append(f"BX:Z:{bx}")
        v = rng.random()
        has_bx = any(t.startswith("BX:") for t in tags)
        if v < 0.05: tags.append("VX:i:0")
        elif v < 0.08 and not has_bx: pass     # a missing VX is only left out where there is no barcode: djinn infers validity from the barcode, Python does not
        elif v < 0.09: tags.append("VX:i:2")
        else: tags.append("VX:i:1")
        if rng.random() < 0.3: tags.append("NM:i:3")
        if rng.random() < 0.3: tags.append("RG:Z:grp")
        rng.shuffle(tags)
        lines.append("\t".join([f"read{len(lines)}", str(flag), cname, str(pos), "60", cig, "=" if flag & 1 else "*",
                                str(pos + 50) if flag & 1 else "0", str(tlen), seq, qual] + tags))
open(out, "w").write("\n".join(lines) + "\n")
if len(sys.argv) > 4:
    with pysam.AlignmentFile(out, "r") as i, pysam.AlignmentFile(sys.argv[4], "wb", template=i) as o:
        for r in i: o.write(r)
