#!/usr/bin/env python3
"""
New-species DIVERSITY amplicon design (Pm, Poc, Pow, Pk).

Rationale
---------
Pf and Pv carry validated, population-genomics-derived microhaplotype panels (Module V, 163 amplicons).
P. malariae, P. ovale (curtisi/wallikeri) and P. knowlesi have little-to-no published population data, so
we CANNOT select loci by measured polymorphism from a single reference genome. Instead we place amplicons
over the domains of canonical antigen genes that are under balancing/diversifying selection and are
documented as polymorphic across EVERY Plasmodium studied — a strong homology prior for diversity.

Method (identical geometry to design_multiexon_amplicon.py):
  1. tblastn the Pf antigen protein vs each species genome -> ortholog HSPs (never span introns).
  2. Map the Pf variable-domain codon window -> species genomic window through the HSPs.
  3. primer3 design flanking that window (conserved flanks, variable core in the read).
  4. in-silico-PCR specificity vs all four new-species orthologs + human off-target check at assembly.

HONEST CAVEAT (surfaced in docs): loci are chosen by homology to known-diverse antigen domains, NOT by
measured intra-species polymorphism. These are discovery-grade diversity/COI markers pending wet-lab
heterozygosity confirmation. Repeat-rich domains (e.g. CSP central repeat, MSP1 block-2) are deliberately
avoided in favour of SNP-dense conserved-flanked domains that microhaplotype cleanly.
"""
import os
import sys, subprocess, csv, os
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import primer3, in_silico_pcr as ip
from design_resistance_amplicon import load_contig, GLOBAL, GENOME
from design_multiexon_amplicon import hsps_for, marker_genomic
from locate_orthologs import load_fasta, rc, translate, attr
import shutil
# Bin dir for BLAST+/mafft — resolved from $MALAFORGE_BIN or PATH (no hardcoded machine path).
BB = (os.environ.get("MALAFORGE_BIN")
      or (os.path.dirname(shutil.which("blastn") or shutil.which("makeblastdb")
                          or shutil.which("mafft") or shutil.which("tblastn") or "") or ""))
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)
import gzip

# Pf antigen locus tags + the codon window of the documented VARIABLE domain (Pf protein coords).
# Windows chosen to sit in SNP-dense, conserved-flanked domains and AVOID repeat/indel tracts.
# Only clean SINGLE-COPY orthologs are used (EBA/DBP/NBP/MSP7 multigene families have no 1:1 Pf
# ortholog and are deliberately excluded to keep tblastn anchoring paralog-safe).
ANTIGENS = {
    "ama1": ("PF3D7_1133400", (200, 284)),   # AMA1 domain I polymorphic cluster (c1L loop), tightened
    "trap": ("PF3D7_1335900", (95, 175)),    # TRAP/SSP2 vWA adhesive core (SNP-dense sub-window)
    "csp":  ("PF3D7_0304600", (302, 372)),   # CSP C-terminal Th2R/Th3R (avoid central repeat)
    "msp1": ("PF3D7_0930300", (1632, 1706)), # MSP1-19 C-terminal double-EGF (SNP-dense; avoid block-2)
}
# Evidence basis per (gene) for the docs — measured diversity in the target species where known.
EVIDENCE = {
    "ama1": "AMA1 Domain I: Pm highest-pi diversifying antigen (PMID 29321885); balancing in Pf/Pv. "
            "Pk EXCLUDED (PkAMA1 purifying, low diversity, PLoS One 10.1371/journal.pone.0124400).",
    "trap": "TRAP/SSP2 vWA: Pm diversifying, 17 aa variants (PMID 29321885); Pk 74 SNPs/29 haplotypes "
            "(10.1186/s12936-018-2423-1). Ovale = homology candidate (unmeasured).",
    "csp":  "CSP C-terminal Th2R/Th3R: Pm 6/9 haplotypes, Asia-Africa Fst 0.40 (PMID 30012172); "
            "ovale dimorphic + Th2R/Th3R (PMID 31170213); Pk pi 0.029 (PMID 26379157).",
    "msp1": "MSP1 C-terminal p19/p42: Pk MSP1-42 18 haplotypes, balancing in Borneo (PMID 30518419); "
            "Pm block-IV clean SNP window nearby (PMID 36114242). Conserved EGF anchor, alignable.",
}
# Per-(gene,taxon) exclusions where the literature shows the domain is NOT informative.
SKIP = {("ama1", "Pk")}   # PkAMA1 under purifying selection -> weak diversity marker
FLANK = 100
MAXPROD = 400          # platform-agnostic ceiling (tailed <=444 -> robust Illumina 2x250 merge)
TAXA = ["Pm", "Pk", "Poc", "Pow"]


def extract_pf_antigens():
    """Translate the antigen genes from Pf3D7.gff by locus tag (generic mRNA->CDS parse)."""
    fa = load_fasta(f"{BASE}/refs/nuclear/Pf3D7.fna")
    want = {lt for lt, _ in ANTIGENS.values()}
    rid2lt = {}; cds = {}
    with gzip.open(f"{BASE}/refs/anno/Pf3D7.gff.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                par = attr(f[8], "Parent"); rid = attr(f[8], "ID")
                if par and par.startswith("gene-") and par[5:] in want and par[5:] not in rid2lt.values():
                    rid2lt[rid] = par[5:]
            elif f[2] == "CDS":
                par = attr(f[8], "Parent")
                if par in rid2lt:
                    cds.setdefault(rid2lt[par], {"contig": f[0], "strand": f[6], "segs": []})
                    cds[rid2lt[par]]["segs"].append((int(f[3]), int(f[4])))
    lt2gene = {lt: g for g, (lt, _) in ANTIGENS.items()}
    prots = {}
    for lt, g in cds.items():
        segs = sorted(g["segs"], reverse=(g["strand"] == "-"))
        nt = ""
        for s, e in segs:
            piece = fa[g["contig"]][s - 1:e]
            nt += rc(piece) if g["strand"] == "-" else piece
        prots[lt2gene[lt]] = translate(nt).rstrip("*")
    return prots


def design(gene, sp, hsps, contig, ortho_regions, lo, hi):
    seq = load_contig(GENOME[sp], contig)
    mr = marker_genomic(hsps, lo, hi)
    if not mr:
        return None
    ms, me = min(mr), max(mr)
    ts = max(0, ms - FLANK); te = min(len(seq), me + FLANK)
    templ = seq[ts:te]; tgt = [ms - ts, max(1, me - ms)]
    g = dict(GLOBAL)
    g['PRIMER_PRODUCT_SIZE_RANGE'] = [[max(180, me - ms + 60), min(MAXPROD, me - ms + 2 * FLANK)]]
    res = primer3.bindings.design_primers(
        {'SEQUENCE_ID': f"{gene}_{sp}", 'SEQUENCE_TEMPLATE': templ, 'SEQUENCE_TARGET': tgt}, g)
    best = None
    for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
        F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
        amps = {}
        for s2, reg in ortho_regions.items():
            if ip.ispcr_bothstrands(F, R, reg, max_mm=4, clean3=3, min_prod=120, max_prod=1600):
                amps[s2] = 1
        cross = [x for x in amps if x != sp]
        cand = {"F": F, "R": R, "prod": res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                "Ftm": round(res[f'PRIMER_LEFT_{i}_TM'], 1), "Rtm": round(res[f'PRIMER_RIGHT_{i}_TM'], 1),
                "self": sp in amps, "ncross": len(cross),
                "spec": "species-specific" if (not cross and sp in amps) else "cross:" + ",".join(cross)}
        if cand["spec"] == "species-specific":
            return cand
        if cand["self"] and (best is None or cand["ncross"] < best["ncross"]):
            best = cand
    return best


def main():
    os.chdir(BASE)
    prots = extract_pf_antigens()
    print("Pf antigen proteins:", {k: len(v) for k, v in prots.items()})
    rows = []
    for gene, (lt, (lo, hi)) in ANTIGENS.items():
        prot = prots.get(gene)
        if not prot:
            print(f"  !! {gene}: Pf protein not extracted (locus {lt})"); continue
        # locate ortholog + build specificity regions per taxon
        hsp_by_sp = {}; ortho_regions = {}
        for sp in TAXA:
            hsps, contig = hsps_for(prot, sp)
            if not hsps:
                print(f"  {sp} {gene}: no ortholog"); continue
            hsp_by_sp[sp] = (hsps, contig)
            ss = [int(h[2]) for h in hsps] + [int(h[3]) for h in hsps]
            ortho_regions[sp] = load_contig(GENOME[sp], contig)[max(0, min(ss) - 300):max(ss) + 300]
        for sp in TAXA:
            if sp not in hsp_by_sp or (gene, sp) in SKIP:
                if (gene, sp) in SKIP:
                    print(f"  {sp} {gene}: SKIP (purifying/low-diversity per literature)")
                continue
            hsps, contig = hsp_by_sp[sp]
            d = design(gene, sp, hsps, contig, ortho_regions, lo, hi)
            aid = f"div_{sp}_{gene}_{lo}-{hi}"
            if d:
                rows.append((sp, gene, aid, d["F"], d["R"], d["prod"], d["Ftm"], d["Rtm"],
                             d["spec"], "homology-anchored", EVIDENCE[gene]))
                print(f"  {sp} {gene} {lo}-{hi}: {d['prod']}bp Tm {d['Ftm']}/{d['Rtm']} [{d['spec']}]")
            else:
                print(f"  {sp} {gene} {lo}-{hi}: NO DESIGN")
    out = f"{BASE}/design/module_V_newsp_diversity.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "gene", "amplicon", "fwd_primer", "rev_primer",
                    "product_bp", "F_tm", "R_tm", "specificity", "grade", "evidence"])
        w.writerows(rows)
    from collections import Counter
    print(f"\n[diversity] {len(rows)} new-species diversity amplicons -> {out}")
    print("  by species:", dict(Counter(r[0] for r in rows)))
    print("  by gene:   ", dict(Counter(r[1] for r in rows)))
    print("  specificity:", dict(Counter('spec' if r[8] == 'species-specific' else 'cross' for r in rows)))


if __name__ == "__main__":
    main()
