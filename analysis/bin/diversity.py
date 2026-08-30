#!/usr/bin/env python3
"""
diversity.py — per-species parasite diversity from microhaplotype amplicons:
  * per-locus number of microhaplotypes (effective alleles) and expected heterozygosity He = 1-Σpᵢ²
  * within-host complexity of infection (COI) — a robust per-locus-count estimator (90th percentile)
  * summarised per species (Pf, Pv, and the new-species antigen loci)

Usage: diversity.py <bam> <manifest.tsv> <sample> <outdir> --species Pf,Pv \
         [--min-depth 25] [--var-maf 0.10] [--hap-freq 0.05]
Outputs: <sample>.diversity.tsv (per-species summary) and <sample>.microhap.tsv (per-locus).
"""
import sys, os, csv, argparse, statistics
from collections import Counter
import numpy as np
import pysam
from common import load_manifest, haplotype_reads


def call_locus(reads, var_maf, hap_freq, min_reads=5):
    """return (n_alleles, He, n_variable_positions, {allele: read_fraction}) from aligned read strings."""
    if len(reads) < min_reads:
        return None
    L = len(reads[0])
    var = []
    for j in range(L):
        col = [r[j] for r in reads if r[j] != "-"]
        if len(col) < min_reads:
            continue
        cnt = Counter(col); tot = sum(cnt.values())
        second = cnt.most_common(2)
        if len(second) > 1 and second[1][1] / tot >= var_maf:
            var.append(j)
    if not var:
        return (1, 0.0, 0, {"REF": 1.0})                  # monomorphic
    haps = Counter("".join(r[j] for j in var) for r in reads if all(r[j] != "-" for j in var))
    tot = sum(haps.values())
    kept = {h: c for h, c in haps.items() if c / tot >= hap_freq and c >= min_reads}
    if not kept:
        return (1, 0.0, len(var), {"REF": 1.0})
    ktot = sum(kept.values())
    fracs = {h: c / ktot for h, c in kept.items()}
    He = 1 - sum(p * p for p in fracs.values())
    return (len(kept), He, len(var), fracs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam"); ap.add_argument("manifest"); ap.add_argument("sample"); ap.add_argument("outdir")
    ap.add_argument("--species", default="Pf,Pv,Poc,Pow,Pm,Pk")
    ap.add_argument("--min-depth", type=int, default=25)
    ap.add_argument("--var-maf", type=float, default=0.10)
    ap.add_argument("--hap-freq", type=float, default=0.05)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    keep = set(a.species.split(","))
    man = load_manifest(a.manifest)
    bam = pysam.AlignmentFile(a.bam, "rb")

    per_locus = []; allele_rows = []; consensus_rows = []
    by_sp = {sp: [] for sp in keep}
    for m in man:
        if m["role"] != "diversity" or m["species"] not in keep:
            continue
        reads = haplotype_reads(bam, m["pan_contig"], m["start"], m["end"])
        if len(reads) < a.min_depth:
            continue
        res = call_locus(reads, a.var_maf, a.hap_freq)
        if not res:
            continue
        nall, He, nvar, fracs = res
        per_locus.append([a.sample, m["species"], m["amplicon"], len(reads), nall, round(He, 3), nvar])
        by_sp[m["species"]].append((nall, He))
        for hap, fr in sorted(fracs.items(), key=lambda x: -x[1]):
            allele_rows.append([a.sample, m["species"], m["amplicon"], hap, round(fr, 3), len(reads)])
        # per-locus majority consensus (population-informative sites are found across samples in relatedness.py)
        L = len(reads[0])
        cons = "".join((Counter(r[j] for r in reads if r[j] != "-").most_common(1) or [("N", 0)])[0][0]
                       for j in range(L))
        consensus_rows.append([a.sample, m["species"], m["amplicon"], cons])

    with open(os.path.join(a.outdir, f"{a.sample}.microhap.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "locus", "depth", "n_microhaplotypes", "He", "n_variable_pos"])
        w.writerows(per_locus)

    # raw per-locus allele calls (within-sample microhaplotypes + read fraction) — input for COI-MLE
    with open(os.path.join(a.outdir, f"{a.sample}.alleles.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "locus", "allele", "read_fraction", "depth"])
        w.writerows(allele_rows)

    # per-locus consensus sequence — cross-sample relatedness uses population-variable positions of these
    with open(os.path.join(a.outdir, f"{a.sample}.consensus.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "locus", "consensus"])
        w.writerows(consensus_rows)

    with open(os.path.join(a.outdir, f"{a.sample}.diversity.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "loci_typed", "polymorphic_loci", "mean_alleles_per_locus",
                    "mean_He", "COI_estimate"])
        for sp in keep:
            vals = by_sp[sp]
            if not vals:
                continue
            counts = [v[0] for v in vals]; hes = [v[1] for v in vals]
            poly = sum(1 for c in counts if c > 1)
            coi = int(max(1, round(np.percentile(counts, 90)))) if counts else 1
            w.writerow([a.sample, sp, len(vals), poly, round(statistics.mean(counts), 2),
                        round(statistics.mean(hes), 3), coi])
            sys.stderr.write(f"[diversity] {a.sample} {sp}: {len(vals)} loci, {poly} polymorphic, COI~{coi}\n")


if __name__ == "__main__":
    main()
