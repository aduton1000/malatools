#!/usr/bin/env python3
"""
generate_barcodes.py — design an orthogonal set of 24-nt sample barcodes for the two-step tailed-amplicon
workflow (round-2 barcoding PCR). Barcodes are chosen for: GC 40-60%, no homopolymer > 3, and pairwise
Hamming distance >= MIN_DIST (default 10) so 2-mismatch demultiplexing is unambiguous. Deterministic.

Usage: generate_barcodes.py [N=96] [LEN=24] [MIN_DIST=10] > barcodes.tsv
Output: TSV  barcode_id<TAB>sequence   (BC001..BCNNN)
"""
import sys, random


def homopolymer_ok(s, maxrun=3):
    run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        if run > maxrun:
            return False
    return True


def gc_ok(s, lo=0.40, hi=0.60):
    gc = (s.count("G") + s.count("C")) / len(s)
    return lo <= gc <= hi


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 96
    L = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    rng = random.Random(20240830)
    chosen = []
    tries = 0
    while len(chosen) < N and tries < 2_000_000:
        tries += 1
        s = "".join(rng.choice("ACGT") for _ in range(L))
        if not gc_ok(s) or not homopolymer_ok(s):
            continue
        if all(hamming(s, c) >= MIN for c in chosen):
            chosen.append(s)
    if len(chosen) < N:
        sys.stderr.write(f"[generate_barcodes] only found {len(chosen)}/{N} at min-dist {MIN}\n")
    sys.stdout.write("barcode_id\tsequence\n")
    for i, s in enumerate(chosen, 1):
        sys.stdout.write(f"BC{i:03d}\t{s}\n")
    sys.stderr.write(f"[generate_barcodes] {len(chosen)} barcodes, {L} nt, min pairwise Hamming {MIN}, "
                     f"GC 40-60%, no homopolymer>3 ({tries} draws)\n")


if __name__ == "__main__":
    main()
