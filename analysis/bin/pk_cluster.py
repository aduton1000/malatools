#!/usr/bin/env python3
"""
pk_cluster.py — assign a P. knowlesi sample to its subpopulation (Mf / Mn / cluster-3) by reading the
fixed cluster-diagnostic SNP alleles in the chr1 C1A, chr6 C2J and chr8 cap380 amplicons, and report the
within-cluster microsatellite alleles.

Diagnostic panel (--panel): amplicon<TAB>offset<TAB>Mf<TAB>Mn<TAB>cluster3  (offset = 0-based position
within the amplicon). In production this is populated from the Divis 2020 fixed-difference SNP list
(10.3201/eid2607.190924) + the chr8 introgression; here it is auto-generated for validation.

Usage: pk_cluster.py <bam> <manifest.tsv> <panel.tsv> <sample> <out.tsv> [--min-depth 15]
"""
import sys, csv, argparse
from collections import Counter, defaultdict
import pysam
from common import load_manifest


def base_at(bam, contig, pos0, min_mapq=10):
    try:
        cov = bam.count_coverage(contig, pos0, pos0 + 1, quality_threshold=0,
                                 read_callback=lambda r: (not r.is_unmapped) and r.mapping_quality >= min_mapq)
    except (ValueError, KeyError):
        return None, 0
    counts = {"A": cov[0][0], "C": cov[1][0], "G": cov[2][0], "T": cov[3][0]}
    d = sum(counts.values())
    if d == 0:
        return None, 0
    return max(counts, key=counts.get), d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam"); ap.add_argument("manifest"); ap.add_argument("panel")
    ap.add_argument("sample"); ap.add_argument("out"); ap.add_argument("--min-depth", type=int, default=15)
    a = ap.parse_args()
    man = {m["amplicon"]: m for m in load_manifest(a.manifest)}
    bam = pysam.AlignmentFile(a.bam, "rb")
    votes = Counter(); loci = []
    for r in csv.DictReader(open(a.panel), delimiter="\t"):
        amp = r["amplicon"]
        if amp not in man:
            continue
        m = man[amp]; pos0 = m["start"] + int(r["offset"])
        base, dep = base_at(bam, m["pan_contig"], pos0)
        if base is None or dep < a.min_depth:
            loci.append([amp, "-", dep, "no_call"]); continue
        callmap = {r["Mf"]: "Mf", r["Mn"]: "Mn", r["cluster3"]: "cluster3"}
        cl = callmap.get(base, "unassigned")
        if cl in ("Mf", "Mn", "cluster3"):
            votes[cl] += 1
        loci.append([amp, base, dep, cl])
    total = sum(votes.values())
    if total == 0:
        call = "no_call"; conf = 0.0
    else:
        call, n = votes.most_common(1)[0]; conf = n / total
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "level", "amplicon", "allele", "depth", "assignment"])
        w.writerow([a.sample, "CALL", "-", "-", "-", f"{call} (confidence {conf:.2f}; votes {dict(votes)})"])
        for amp, base, dep, cl in loci:
            w.writerow([a.sample, "locus", amp, base, dep, cl])
    sys.stderr.write(f"[pk_cluster] {a.sample}: cluster = {call} (conf {conf:.2f}, votes {dict(votes)})\n")


if __name__ == "__main__":
    main()
