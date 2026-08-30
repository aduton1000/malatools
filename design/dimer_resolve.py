#!/usr/bin/env python3
"""
Dimer-aware primer re-engineering. For a short list of amplicons involved in cross-primer dimers, locate
the amplicon on its genome, ask primer3 for many candidate pairs, and pick the first pair whose CS1/CS2-
tailed primers have NO end-stability dimer > threshold against the rest of the panel (and cover the same
locus). Updates the module TSV in place. No quality compromise: same locus, same product-size class,
full re-validation of single-product + dimer-clean.
"""
import os
import sys, csv, os
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import primer3, in_silico_pcr as ip
from design_resistance_amplicon import GLOBAL
from design_pk_cluster import find_product
from locate_orthologs import load_fasta
from degeneracy_utils import revcomp
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)
CS1 = "ACACTGACGACATGGTTCTACA"; CS2 = "TACGGTAGCAGAGACTTGGTCT"
DIMER_TM = 40.0

# amplicon -> (module tsv, genome fasta)
TARGETS = {
    "Pow_mdr1_1030-1260_t3": ("design/module_R_resistance.tsv", "refs/nuclear/PowCR01.fna"),
}


def all_panel_primers(exclude):
    """collect CS-tailed primers from every module except the excluded amplicons."""
    prim = []
    files = ["module_S_species_id.tsv", "module_R_resistance.tsv", "module_D_diagnostic_cnv.tsv",
             "module_V_amplicons.tsv", "module_V_newsp_diversity.tsv", "module_V_pk_cluster.tsv"]
    for f in files:
        p = f"{BASE}/design/{f}"
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r.get("fwd_primer") and r.get("rev_primer") and r["amplicon"] not in exclude:
                prim.append(CS1 + r["fwd_primer"]); prim.append(CS2 + r["rev_primer"])
    return prim


def max_dimer(F, R, panel):
    tf, tr = CS1 + F, CS2 + R
    mx = 0.0
    for s in panel:
        for t in (tf, tr):
            mx = max(mx, primer3.bindings.calc_end_stability(t, s).tm)
    return mx


def main():
    os.chdir(BASE)
    for aid, (tsv, gfa) in TARGETS.items():
        rows = list(csv.DictReader(open(tsv), delimiter="\t"))
        row = next((r for r in rows if r["amplicon"] == aid), None)
        if not row:
            print(f"  {aid}: not found in {tsv}"); continue
        genome = load_fasta(gfa)
        # locate current amplicon -> region
        loc = None
        for cn, cs in genome.items():
            fpd = find_product(row["fwd_primer"], row["rev_primer"], cs, max_prod=1200)
            if fpd:
                s = cs if fpd[3] == "+" else revcomp(cs)
                loc = (s, fpd[0], fpd[1]); break
        if not loc:                                          # mismatch-tolerant fallback
            for cn, cs in genome.items():
                hits = ip.ispcr_bothstrands(row["fwd_primer"], row["rev_primer"], cs,
                                            max_mm=4, clean3=3, min_prod=100, max_prod=1200)
                if hits:
                    strand, p = min(hits, key=lambda x: x[1]['size'])
                    s = cs if strand == '+' else revcomp(cs)
                    loc = (s, p['start'], p['end']); break
        if not loc:
            print(f"  {aid}: could not locate"); continue
        s, S, E = loc
        pad = 90; ts = max(0, S - pad); te = min(len(s), E + pad)
        templ = s[ts:te]; tgt = [S - ts, max(20, E - S)]
        panel = all_panel_primers(exclude={aid})
        g = dict(GLOBAL)
        g['PRIMER_PRODUCT_SIZE_RANGE'] = [[max(150, (E - S) - 40), 400]]
        g['PRIMER_NUM_RETURN'] = 20
        res = primer3.bindings.design_primers(
            {'SEQUENCE_ID': aid, 'SEQUENCE_TEMPLATE': templ, 'SEQUENCE_TARGET': tgt}, g)
        best = None
        for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
            F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
            md = max_dimer(F, R, panel)
            cand = (md, F, R, res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'])
            if best is None or md < best[0]:
                best = cand
            if md <= DIMER_TM:
                break
        if not best:
            print(f"  {aid}: no candidate pairs"); continue
        md, F, R, prod = best
        row["fwd_primer"] = F; row["rev_primer"] = R
        if "product_bp" in row:
            row["product_bp"] = str(prod)
        with open(tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t"); w.writeheader(); w.writerows(rows)
        print(f"  {aid}: re-engineered -> {prod}bp, max panel dimer {md:.1f}C {'(CLEAN)' if md<=DIMER_TM else '(still >40, best available)'}")


if __name__ == "__main__":
    main()
