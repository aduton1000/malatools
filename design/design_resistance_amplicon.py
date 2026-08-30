#!/usr/bin/env python3
"""
Design species-specific resistance amplicons over a gene's marker-codon region, and validate:
- primer3 design on the species genome (AT-rich-tuned)
- in-silico PCR of the chosen pair against EVERY species' ortholog region -> species-specificity
- confirm the marker region is inside the insert

Marker windows (species codon numbering, generous to cover orthologs):
  dhfr: 8-135   (antifolate markers 15/51/57/58/108/113/114 ...)
  k13 : 335-560 (propeller; artemisinin discovery)
Usage: design_resistance_amplicon.py <gene>
"""
import os
import sys, gzip, subprocess
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import primer3, in_silico_pcr as ip
from degeneracy_utils import revcomp
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)
GENOME = {"Pm": "PmUG01", "Pk": "PkHl", "Poc": "PocGH01", "Pow": "PowCR01", "Pv": "PvP01", "Pf": "Pf3D7"}
MARKER = {"dhfr": (12, 120), "k13": (335, 560), "dhps": (390, 620), "crt": (60, 380), "mdr1": (70, 260)}
FLANK = 230
GLOBAL = {'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 28, 'PRIMER_OPT_TM': 59,
          'PRIMER_MIN_TM': 52, 'PRIMER_MAX_TM': 64, 'PRIMER_PAIR_MAX_DIFF_TM': 6, 'PRIMER_MIN_GC': 15,
          'PRIMER_OPT_GC_PERCENT': 30, 'PRIMER_MAX_GC': 70, 'PRIMER_MAX_POLY_X': 6,
          'PRIMER_MAX_NS_ACCEPTED': 0, 'PRIMER_NUM_RETURN': 12, 'PRIMER_MAX_SELF_ANY_TH': 47,
          'PRIMER_MAX_SELF_END_TH': 40, 'PRIMER_MAX_HAIRPIN_TH': 47}

def load_contig(genome, contig):
    d = {}; name = None; buf = []
    for line in open(f"{BASE}/refs/nuclear/{genome}.fna"):
        if line.startswith(">"):
            if name == contig: return "".join(buf).upper()
            name = line[1:].split()[0]; buf = []
        elif name == contig:
            buf.append(line.strip())
    return "".join(buf).upper()

def orthologs():
    rows = [l.split("\t") for l in open(f"{BASE}/targets/orthologs.tsv").read().splitlines()[1:]]
    return {(r[0], r[1]): (r[2], int(r[3]), int(r[4]), r[5]) for r in rows}

def marker_region(gene, start, end, strand):
    lo, hi = MARKER[gene]
    if strand == "+":
        return start + 3 * (lo - 1), start + 3 * hi
    else:
        return end - 3 * hi, end - 3 * (lo - 1)

def design_for(gene, sp, orth, ortho_regions):
    contig, start, end, strand = orth[(gene, sp)]
    seq = load_contig(GENOME[sp], contig)
    ms, me = marker_region(gene, start, end, strand)
    ts = max(0, ms - FLANK); te = min(len(seq), me + FLANK)
    templ = seq[ts:te]
    tgt = [ms - ts, me - ms]
    g = dict(GLOBAL); g['PRIMER_PRODUCT_SIZE_RANGE'] = [[me - ms + 60, me - ms + 2 * FLANK]]
    res = primer3.bindings.design_primers({'SEQUENCE_ID': f"{gene}_{sp}", 'SEQUENCE_TEMPLATE': templ,
                                           'SEQUENCE_TARGET': tgt}, g)
    best = None
    for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
        F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
        # species-specificity: in-silico PCR against every ortholog region
        amps = {}
        for s2, reg in ortho_regions.items():
            hits = ip.ispcr_bothstrands(F, R, reg, max_mm=4, clean3=3, min_prod=120, max_prod=1500)
            if hits: amps[s2] = min(h[1]["size"] for h in hits)
        cross = [x for x in amps if x != sp]
        cand = {"sp": sp, "gene": gene, "F": F, "R": R, "prod": res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                "Ftm": round(res[f'PRIMER_LEFT_{i}_TM'], 1), "Rtm": round(res[f'PRIMER_RIGHT_{i}_TM'], 1),
                "ncross": len(cross), "amps_self": sp in amps,
                "spec": "species-specific" if not cross and sp in amps else "cross:" + ",".join(cross)}
        if cand["spec"] == "species-specific":
            return cand
        # keep the most-specific pair that still amplifies the target species
        if cand["amps_self"] and (best is None or cand["ncross"] < best["ncross"]):
            best = cand
    return best

def main():
    gene = sys.argv[1]
    orth = orthologs()
    # extract each species' ortholog region (gene +/- 300 bp) for specificity testing
    ortho_regions = {}
    for (g, sp), (contig, start, end, strand) in orth.items():
        if g != gene: continue
        seq = load_contig(GENOME[sp], contig)
        ortho_regions[sp] = seq[max(0, start - 300):end + 300]
    print(f"gene={gene}  species with ortholog: {list(ortho_regions)}")
    out = []
    for sp in ["Pm", "Pk", "Poc", "Pow"]:
        if (gene, sp) not in orth: continue
        r = design_for(gene, sp, orth, ortho_regions)
        if r:
            out.append(r)
            print(f"  {sp} {gene}: prod {r['prod']}bp  Tm {r['Ftm']}/{r['Rtm']}  [{r['spec']}]")
            print(f"       F {r['F']}\n       R {r['R']}")
        else:
            print(f"  {sp} {gene}: NO DESIGN")
    with open(f"{BASE}/design/resistance_{gene}.tsv", "w") as fh:
        fh.write("species\tgene\tfwd_primer\trev_primer\tproduct_bp\tF_tm\tR_tm\tspecificity\n")
        for r in out:
            fh.write(f"{r['sp']}\t{r['gene']}\t{r['F']}\t{r['R']}\t{r['prod']}\t{r['Ftm']}\t{r['Rtm']}\t{r['spec']}\n")

if __name__ == "__main__":
    main()
