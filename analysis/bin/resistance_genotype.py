#!/usr/bin/env python3
"""
resistance_genotype.py — call drug-resistance genotypes per species at every marker codon, ALWAYS vs the
SENSITIVE (wild-type) allele so reference-is-mutant loci (e.g. pmdhfr 58/114) are not silently miscalled.
Also reports discovery-grade non-synonymous variants in the new-species resistance amplicons.

Usage: resistance_genotype.py <bam> <manifest.tsv> <codon_map.tsv> <sample> <outdir> \
         --species Pf,Pv [--min-depth 10] [--min-vaf 0.05]
Outputs: <outdir>/<sample>.resistance.tsv (per-codon calls) and <sample>.drugs.tsv (per-drug summary).
"""
import sys, os, csv, argparse
from collections import defaultdict
import pysam
from common import load_manifest, codon_alleles, species_of, pan_contig


def load_codons(path):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        try:
            r["codon"] = int(r["codon"]); r["genomic_start"] = int(r["genomic_start"]); r["genomic_end"] = int(r["genomic_end"])
        except (ValueError, TypeError):
            continue
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam"); ap.add_argument("manifest"); ap.add_argument("codon_map")
    ap.add_argument("sample"); ap.add_argument("outdir")
    ap.add_argument("--species", default="Pf,Pv,Poc,Pow,Pm,Pk")
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-vaf", type=float, default=0.05)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    keep_sp = set(a.species.split(","))
    codons = load_codons(a.codon_map)
    bam = pysam.AlignmentFile(a.bam, "rb")

    rows = []
    drug_hits = defaultdict(lambda: defaultdict(list))    # (sp) -> drug -> [mutant marker strings]
    for c in codons:
        sp = species_of(c["gene"])
        if sp not in keep_sp:
            continue
        pc = pan_contig(sp, c["contig"])
        alleles, depth = codon_alleles(bam, pc, c["genomic_start"], c["genomic_end"], c["strand"])
        wt = (c.get("wt_aa") or "").strip()
        drugs = (c.get("drugs") or "").strip()
        if depth < a.min_depth:
            call = "no_call(depth<%d)" % a.min_depth
            major = "-"; alts = ""
        else:
            major = max(alleles, key=alleles.get)
            if wt:
                muts = {aa: f for aa, f in alleles.items() if aa != wt and f >= a.min_vaf}
                call = "RESISTANT" if muts else "sensitive"
                alts = ";".join(f"{aa}:{f:.2f}" for aa, f in sorted(muts.items(), key=lambda x: -x[1]))
                if muts:
                    for d in drugs.split("/"):
                        drug_hits[sp][d.strip()].append(f"{c['gene']}:{wt}{c['codon']}{major}")
            else:
                call = "genotype"; alts = f"{major}:{alleles[major]:.2f}"
        rows.append([a.sample, sp, c["gene"], c["codon"], c.get("ref_aa", ""), wt or "?", major,
                     alts, depth, call, drugs])

    rtsv = os.path.join(a.outdir, f"{a.sample}.resistance.tsv")
    with open(rtsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "gene", "codon", "ref_aa", "sensitive_aa", "observed_aa",
                    "resistant_alleles(vaf)", "depth", "call", "drugs"])
        w.writerows(rows)

    dtsv = os.path.join(a.outdir, f"{a.sample}.drugs.tsv")
    with open(dtsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "drug", "status", "markers"])
        drugset = set()
        for sp in keep_sp:
            for c in codons:
                if species_of(c["gene"]) == sp:
                    for d in (c.get("drugs") or "").split("/"):
                        if d.strip():
                            drugset.add((sp, d.strip()))
        for sp, d in sorted(drugset):
            hits = drug_hits[sp].get(d, [])
            status = "RESISTANCE-ASSOCIATED MUTATION" if hits else "no resistance marker detected"
            w.writerow([a.sample, sp, d, status, ";".join(sorted(set(hits)))])
    nres = sum(1 for r in rows if r[9] == "RESISTANT")
    sys.stderr.write(f"[resistance] {a.sample}: {nres} resistant codon call(s) across {len(keep_sp)} species -> {rtsv}\n")


if __name__ == "__main__":
    main()
