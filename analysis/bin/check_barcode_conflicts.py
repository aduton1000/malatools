#!/usr/bin/env python3
"""
check_barcode_conflicts.py — verify the sample-barcode set is consistent with the rest of the library
construct and the target genomes. Run this whenever the barcode set changes (e.g. swapping in reagent
barcodes). Exits non-zero if any check fails.

Checks:
  1. barcode set is internally orthogonal (min pairwise Hamming >= --min-dist)
  2. barcodes are dissimilar from CS1 / CS2 (so CS-anchored demux is unambiguous)
  3. no barcode is within --primer-mm of any panel core-primer window (mispriming risk)
  4. CS1 / CS2 do not occur in the pan-reference (CS tails never prime on target DNA)

Usage: check_barcode_conflicts.py [--barcodes refs/barcodes/barcodes_96.tsv] [--min-dist 8] [--primer-mm 3]
"""
import os, sys, csv, gzip, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.abspath(os.path.join(HERE, "../.."))
CS1 = "ACACTGACGACATGGTTCTACA"; CS2 = "TACGGTAGCAGAGACTTGGTCT"
COMP = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(COMP)[::-1]
def ham(a, b): return sum(1 for x, y in zip(a, b) if x != y)
def minwin(short, long): return min((ham(short, long[i:i + len(short)]) for i in range(len(long) - len(short) + 1)), default=99)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcodes", default=f"{BASE}/refs/barcodes/barcodes_96.tsv")
    ap.add_argument("--min-dist", type=int, default=8)
    ap.add_argument("--primer-mm", type=int, default=3)
    ap.add_argument("--ref", default=f"{BASE}/refs/pan_reference.fa.gz")
    a = ap.parse_args()
    bc = {r["barcode_id"]: r["sequence"].upper() for r in csv.DictReader(open(a.barcodes), delimiter="\t")}
    ids = list(bc); seqs = [bc[i] for i in ids]
    fails = []

    # 1. orthogonality
    mn = min(ham(seqs[i], seqs[j]) for i in range(len(seqs)) for j in range(i + 1, len(seqs)))
    ok1 = mn >= 10
    print(f"1. barcode orthogonality: min pairwise Hamming = {mn}  [{'OK' if ok1 else 'FAIL'}]")
    if not ok1: fails.append("orthogonality")

    # 2. vs CS1/CS2
    worst = min(min(minwin(CS1, s + s), minwin(CS2, s + s), minwin(CS1, rc(s) + rc(s)), minwin(CS2, rc(s) + rc(s)))
                for s in seqs)
    ok2 = worst >= a.min_dist
    print(f"2. barcodes vs CS1/CS2: worst-case window Hamming = {worst} (>= {a.min_dist})  [{'OK' if ok2 else 'FAIL'}]")
    if not ok2: fails.append("vs-CS")

    # 3. vs panel primers
    prims = set()
    for f in ["module_S_species_id", "module_R_resistance", "module_D_diagnostic_cnv",
              "module_V_amplicons", "module_V_newsp_diversity", "module_V_pk_cluster"]:
        p = f"{BASE}/design/{f}.tsv"
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                for k in ("fwd_primer", "rev_primer"):
                    if r.get(k): prims.add(r[k].upper())
    hits = []
    for bid, s in bc.items():
        for p in prims:
            long, short = (p, s) if len(p) >= len(s) else (s, p)
            if minwin(short, long) <= a.primer_mm:
                hits.append((bid, p)); break
    ok3 = not hits
    print(f"3. barcodes vs {len(prims)} panel primers: {len(hits)} within {a.primer_mm} mm  [{'OK' if ok3 else 'FAIL'}]")
    if not ok3: fails.append("vs-primers")

    # 4. CS1/CS2 in reference
    tags = {CS1: 0, CS2: 0, rc(CS1): 0, rc(CS2): 0}
    buf = ""
    with gzip.open(a.ref, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                for t in tags: tags[t] += buf.count(t)
                buf = ""
            else:
                buf += line.strip()
                if len(buf) > 200000:
                    for t in tags: tags[t] += buf[:-30].count(t)
                    buf = buf[-30:]
        for t in tags: tags[t] += buf.count(t)
    ncs = sum(tags.values())
    ok4 = ncs == 0
    print(f"4. CS1/CS2 in pan-reference: {ncs} occurrences (all orientations)  [{'OK' if ok4 else 'FAIL'}]")
    if not ok4: fails.append("CS-in-genome")

    print("\n" + ("✓ ALL CHECKS PASS — barcodes are consistent with the library construct" if not fails
                  else "✗ FAILED: " + ", ".join(fails)))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
