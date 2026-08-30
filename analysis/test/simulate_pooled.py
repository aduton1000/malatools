#!/usr/bin/env python3
"""
simulate_pooled.py — build a POOLED, barcoded ONT-style FASTQ from the per-sample amplicon reads, to test
demux.py. Each read is wrapped:  [adapter][fwd_bc][CS1][amplicon][revcomp CS2][revcomp rev_bc][adapter],
randomly reverse-complemented (ONT reads both strands), with substitution errors. Also writes the
sample->barcode map that demux.py consumes.

Usage: simulate_pooled.py <reads_dir> <barcodes.tsv> <pooled.fastq.gz> <map.csv>
"""
import sys, os, gzip, csv, glob, random

CS1 = "ACACTGACGACATGGTTCTACA"; CS2 = "TACGGTAGCAGAGACTTGGTCT"
COMP = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(COMP)[::-1]
rng = random.Random(7)
ERR = 0.02          # ONT-ish
ADAP = "AATGTACTTCGTTCAGTTACGTATTGCT"   # a stand-in leader/adapter fragment


def mutate(s):
    b = list(s)
    for i in range(len(b)):
        if rng.random() < ERR:
            b[i] = rng.choice([x for x in "ACGT" if x != b[i]])
    return "".join(b)


def main():
    reads_dir, bc_tsv, pooled, map_out = sys.argv[1:5]
    bc = {r["barcode_id"]: r["sequence"] for r in csv.DictReader(open(bc_tsv), delimiter="\t")}
    ids = list(bc)
    samples = sorted(os.path.basename(f)[:-9] for f in glob.glob(f"{reads_dir}/*.fastq.gz"))
    # assign each sample a distinct (fwd,rev) barcode pair
    assign = {}
    for i, s in enumerate(samples):
        assign[s] = (ids[i % len(ids)], ids[(i + len(samples)) % len(ids)])
    with open(map_out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["sample", "fwd_bc", "rev_bc"])
        for s in samples:
            w.writerow([s, assign[s][0], assign[s][1]])

    n = 0
    with gzip.open(pooled, "wt") as out:
        for s in samples:
            fbid, rbid = assign[s]; fseq, rseq = bc[fbid], bc[rbid]
            with gzip.open(f"{reads_dir}/{s}.fastq.gz", "rt") as fh:
                while True:
                    h = fh.readline()
                    if not h:
                        break
                    amp = fh.readline().strip(); fh.readline(); fh.readline()
                    mol = ADAP[rng.randint(0, 12):] + fseq + CS1 + amp + rc(CS2) + rc(rseq) + ADAP[:rng.randint(6, 20)]
                    mol = mutate(mol)
                    if rng.random() < 0.5:
                        mol = rc(mol)
                    out.write(f"{h.strip()}\n{mol}\n+\n{'I'*len(mol)}\n"); n += 1
    print(f"[simulate_pooled] {n} barcoded reads across {len(samples)} samples -> {pooled}")
    print(f"[simulate_pooled] map -> {map_out}")


if __name__ == "__main__":
    main()
