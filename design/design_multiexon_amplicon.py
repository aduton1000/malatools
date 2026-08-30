#!/usr/bin/env python3
"""
Exon-aware resistance-amplicon design (for crt, mdr1, and any multi-exon gene).

Maps Pf marker codons to each species' genome through tblastn HSPs (which never span introns), so the
amplicon is placed on the true marker-bearing exon(s). Then primer3 + in-silico-PCR specificity, as in
design_resistance_amplicon.py.

Target Pf-codon windows (WHO/known marker regions):
  crt : [72,150]                     (K76T / CVIET-SVMNT region + neighbours)
  mdr1: [80,200] and [1030,1260]     (N86/Y184 ; and S1034/N1042/D1246 region) -> two amplicons
"""
import os
import sys, subprocess
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import primer3, in_silico_pcr as ip
from design_resistance_amplicon import load_contig, GLOBAL, GENOME, orthologs
from locate_orthologs import extract_pf_proteins
import shutil
# Bin dir for BLAST+/mafft — resolved from $MALAFORGE_BIN or PATH (no hardcoded machine path).
BB = (os.environ.get("MALAFORGE_BIN")
      or (os.path.dirname(shutil.which("blastn") or shutil.which("makeblastdb")
                          or shutil.which("mafft") or shutil.which("tblastn") or "") or ""))
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)
WINDOWS = {"crt": [(72, 150)], "mdr1": [(80, 200), (1030, 1260)]}
FLANK = 220

def hsps_for(prot, sp):
    open("scratch/q.faa", "w").write(f">q\n{prot}\n")
    r = subprocess.run([f"{BB}/tblastn", "-query", "scratch/q.faa", "-db", f"scratch/db_{sp}",
                        "-outfmt", "6 qstart qend sstart send sframe sseqid pident", "-max_target_seqs", "8"],
                       capture_output=True, text=True)
    hs = [l.split("\t") for l in r.stdout.strip().splitlines()]
    if not hs: return None, None
    # anchor on best HSP; keep same-contig HSPs within 10kb
    tc = hs[0][5]; tmid = (int(hs[0][2]) + int(hs[0][3])) / 2
    keep = [h for h in hs if h[5] == tc and abs((int(h[2]) + int(h[3])) / 2 - tmid) < 12000]
    return keep, tc

def codon_to_genomic(hsps, codon):
    """map a Pf codon to a genomic position via the HSP that contains it (linear within an HSP)."""
    for qs, qe, ss, se, frame, sid, pid in [(int(h[0]), int(h[1]), int(h[2]), int(h[3]), int(h[4]), h[5], h[6]) for h in hsps]:
        if qs <= codon <= qe:
            off = (codon - qs) * 3
            return ss + off if frame > 0 else ss - off
    return None

def marker_genomic(hsps, lo, hi):
    """genomic span covering Pf codons lo..hi (using nearest mappable codons if endpoints fall in introns)."""
    pos = []
    for c in range(lo, hi + 1):
        g = codon_to_genomic(hsps, c)
        if g is not None: pos.append(g)
    if not pos: return None
    return min(pos), max(pos)

def design(gene, sp, hsps, contig, ortho_regions, lo, hi):
    seq = load_contig(GENOME[sp], contig)
    mr = marker_genomic(hsps, lo, hi)
    if not mr: return None
    ms, me = mr; ms, me = min(ms, me), max(ms, me)
    ts = max(0, ms - FLANK); te = min(len(seq), me + FLANK)
    templ = seq[ts:te]; tgt = [ms - ts, max(1, me - ms)]
    g = dict(GLOBAL); g['PRIMER_PRODUCT_SIZE_RANGE'] = [[max(160, me - ms + 60), me - ms + 2 * FLANK]]
    res = primer3.bindings.design_primers({'SEQUENCE_ID': f"{gene}_{sp}", 'SEQUENCE_TEMPLATE': templ,
                                           'SEQUENCE_TARGET': tgt}, g)
    best = None
    for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
        F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
        amps = {}
        for s2, reg in ortho_regions.items():
            if ip.ispcr_bothstrands(F, R, reg, max_mm=4, clean3=3, min_prod=120, max_prod=1600): amps[s2] = 1
        cross = [x for x in amps if x != sp]
        cand = {"F": F, "R": R, "prod": res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                "Ftm": round(res[f'PRIMER_LEFT_{i}_TM'], 1), "Rtm": round(res[f'PRIMER_RIGHT_{i}_TM'], 1),
                "ncross": len(cross), "self": sp in amps,
                "spec": "species-specific" if (not cross and sp in amps) else "cross:" + ",".join(cross)}
        if cand["spec"] == "species-specific": return cand
        if cand["self"] and (best is None or cand["ncross"] < best["ncross"]): best = cand
    return best

def main():
    gene = sys.argv[1]
    prots = extract_pf_proteins(); prot = prots[gene]; orth = orthologs()
    ortho_regions = {}
    for (g, sp), (contig, start, end, strand) in orth.items():
        if g == gene: ortho_regions[sp] = load_contig(GENOME[sp], contig)[max(0, start - 300):end + 300]
    rows = []
    for sp in ["Pm", "Pk", "Poc", "Pow"]:
        hsps, contig = hsps_for(prot, sp)
        if not hsps: print(f"  {sp} {gene}: no ortholog"); continue
        for (lo, hi) in WINDOWS[gene]:
            tag = f"{gene}_{lo}-{hi}"
            d = design(gene, sp, hsps, contig, ortho_regions, lo, hi)
            if d:
                rows.append((sp, tag, d["F"], d["R"], d["prod"], d["Ftm"], d["Rtm"], d["spec"]))
                print(f"  {sp} {tag}: {d['prod']}bp Tm {d['Ftm']}/{d['Rtm']} [{d['spec']}]")
            else:
                print(f"  {sp} {tag}: NO DESIGN")
    with open(f"{BASE}/design/resistance_{gene}.tsv", "w") as fh:
        fh.write("species\tgene\tfwd_primer\trev_primer\tproduct_bp\tF_tm\tR_tm\tspecificity\n")
        for r in rows: fh.write("\t".join(map(str, r)) + "\n")

if __name__ == "__main__":
    main()
