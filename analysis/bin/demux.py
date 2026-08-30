#!/usr/bin/env python3
"""
demux.py — demultiplex a pooled ONT/PCR run by the panel's CUSTOM combinatorial 24-nt barcodes.

MinKNOW's native demux only knows ONT's native barcode kits; this panel uses custom dual 24-nt PCR
barcodes (round-2 barcoding PCR), so demultiplexing is done here. Reads are anchored on the constant CS1
/ CS2 tags (robust to the variable adapter length before the barcode): the 24 nt immediately 5' of CS1
is the forward barcode; the 24 nt immediately 3' of the reverse tag (revcomp CS2) is the reverse barcode.
Both read orientations are tried. Barcodes are matched with mismatch tolerance; a read is assigned to a
sample only when its (fwd,rev) barcode pair is in the sample map. The emitted per-sample FASTQ contains
the trimmed amplicon (CS1…CS2 removed) — ready for MAP.

Usage:
  demux.py --fastq pooled.fastq.gz --barcodes barcodes_96.tsv --sample-map map.csv --outdir demuxed \
     [--mode dual|single] [--bc-mismatch 2] [--cs-mismatch 3] [--search 60] [--min-amplicon 120] [--write-samplesheet]
map.csv columns:  sample,fwd_bc[,rev_bc]   (barcode IDs from barcodes.tsv; rev_bc required for --mode dual)
"""
import sys, os, csv, gzip, argparse
from collections import Counter

CS1 = "ACACTGACGACATGGTTCTACA"
CS2 = "TACGGTAGCAGAGACTTGGTCT"
COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(COMP)[::-1]


def load_barcodes(path):
    bc = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        bc[r["barcode_id"]] = r["sequence"].upper()
    return bc


def find_tag(read, tag, max_mm, start, end):
    """min-mismatch start position of `tag`, searching start positions in [start, end); or None."""
    L = len(tag); best = None; bestmm = max_mm + 1
    hi = min(end, len(read) - L + 1)
    i = max(0, start)
    while i < hi:
        mm = 0
        for k in range(L):
            if read[i + k] != tag[k]:
                mm += 1
                if mm > bestmm:
                    break
        if mm < bestmm:
            bestmm = mm; best = i
            if mm == 0:
                break
        i += 1
    return best


def match_bc(window, exact, seqs, max_mm):
    """return barcode_id for a 24-nt window (exact dict first, else nearest within max_mm & unambiguous)."""
    if window in exact:
        return exact[window]
    best = None; bestmm = max_mm + 1; second = max_mm + 1
    for bid, seq in seqs:
        mm = sum(1 for a, b in zip(window, seq) if a != b)
        if mm < bestmm:
            second = bestmm; bestmm, best = mm, bid
        elif mm < second:
            second = mm
    if bestmm <= max_mm and second - bestmm >= 1:      # unambiguous (>=1 nt clearer than runner-up)
        return best
    return None


def fastq_reader(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip("\n"); fh.readline(); q = fh.readline().rstrip("\n")
            yield h.rstrip("\n"), s.upper(), q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True, nargs="+")
    ap.add_argument("--barcodes", required=True); ap.add_argument("--sample-map", required=True)
    ap.add_argument("--outdir", required=True); ap.add_argument("--mode", default="dual", choices=["dual", "single"])
    ap.add_argument("--bc-mismatch", type=int, default=2); ap.add_argument("--cs-mismatch", type=int, default=3)
    ap.add_argument("--search", type=int, default=60); ap.add_argument("--min-amplicon", type=int, default=120)
    ap.add_argument("--write-samplesheet", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    bc = load_barcodes(a.barcodes)
    exact = {v: k for k, v in bc.items()}
    seqs = list(bc.items())
    # sample map keyed by barcode-id pair
    smap = {}
    for r in csv.DictReader(open(a.sample_map)):
        if a.mode == "dual":
            smap[(r["fwd_bc"].strip(), r["rev_bc"].strip())] = r["sample"].strip()
        else:
            smap[r["fwd_bc"].strip()] = r["sample"].strip()
    RCCS2 = rc(CS2)

    handles = {s: gzip.open(f"{a.outdir}/{s}.fastq.gz", "wt") for s in set(smap.values())}
    unc = gzip.open(f"{a.outdir}/unclassified.fastq.gz", "wt")
    counts = Counter(); reasons = Counter(); total = 0

    for fq in a.fastq:
        for h, seq, qual in fastq_reader(fq):
            total += 1
            assigned = None; amp = ampq = None
            for orient in (0, 1):
                r = seq if orient == 0 else rc(seq)
                q = qual if orient == 0 else qual[::-1]
                p1 = find_tag(r, CS1, a.cs_mismatch, 0, a.search)
                if p1 is None or p1 < 24:
                    continue
                p2 = find_tag(r, RCCS2, a.cs_mismatch, p1 + len(CS1) + a.min_amplicon - 20, len(r))
                if p2 is None:
                    continue
                fbid = match_bc(r[p1 - 24:p1], exact, seqs, a.bc_mismatch)
                rbid = match_bc(rc(r[p2 + len(CS2):p2 + len(CS2) + 24]), exact, seqs, a.bc_mismatch)
                key = (fbid, rbid) if a.mode == "dual" else fbid
                if key in smap:
                    assigned = smap[key]; amp = r[p1 + len(CS1):p2]; ampq = q[p1 + len(CS1):p2]
                    break
                reasons["bc_pair_not_in_map" if (fbid or rbid) else "no_barcode"] += 1
            if assigned and len(amp) >= a.min_amplicon:
                handles[assigned].write(f"{h}\n{amp}\n+\n{ampq}\n"); counts[assigned] += 1
            else:
                unc.write(f"{h}\n{seq}\n+\n{qual}\n"); counts["unclassified"] += 1
    for fh in handles.values():
        fh.close()
    unc.close()

    with open(f"{a.outdir}/demux_stats.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["sample", "reads", "pct_of_total"])
        for s in sorted(set(smap.values())):
            w.writerow([s, counts[s], round(100 * counts[s] / max(1, total), 2)])
        w.writerow(["unclassified", counts["unclassified"], round(100 * counts["unclassified"] / max(1, total), 2)])
    assigned_tot = total - counts["unclassified"]
    sys.stderr.write(f"[demux] {total} reads -> {assigned_tot} assigned ({100*assigned_tot/max(1,total):.1f}%), "
                     f"{counts['unclassified']} unclassified across {len(set(smap.values()))} samples\n")

    if a.write_samplesheet:
        with open(f"{a.outdir}/samples.csv", "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["sample", "fastq"])
            for s in sorted(set(smap.values())):
                if counts[s] > 0:
                    w.writerow([s, os.path.abspath(f"{a.outdir}/{s}.fastq.gz")])
        sys.stderr.write(f"[demux] samplesheet -> {a.outdir}/samples.csv\n")


if __name__ == "__main__":
    main()
