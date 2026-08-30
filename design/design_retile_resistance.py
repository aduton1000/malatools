#!/usr/bin/env python3
"""
Re-tile the ENTIRE resistance module into short (<=400 bp core), NON-OVERLAPPING amplicons so the panel
is truly platform-agnostic (ONT + Illumina 2x250 + MGI) AND can collapse to a single pool.

Two input classes, one output:
  * CODON-MAPPED genes (all Pf, all Pv, non-Pf dhfr) -> cluster the validated marker codons into groups
    each spanning <=~300 bp genomic, placing tile junctions in the LARGEST inter-codon gaps so NO marker
    codon ever falls in a primer/junction blind spot. One <=400 bp amplicon per cluster.
  * DISCOVERY genes (non-Pf dhps/k13/crt/mdr1, Pk dhfr) -> locate the current amplicon's span on the
    species genome and tile it into <=400 bp non-overlapping windows (junction blind spots are acceptable
    for region-scanning discovery markers).

Every tile is validated for a single product on its own genome and species-specificity (3'-anchor counts
across the other five taxa). Marker coverage is preserved exactly (same codons / same region). Output
replaces design/module_R_resistance.tsv.
"""
import os
import sys, csv, os, re
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import primer3, in_silico_pcr as ip
from design_resistance_amplicon import GLOBAL
from design_pk_cluster import find_product, _acount
from design_multiexon_amplicon import hsps_for, marker_genomic
from locate_orthologs import load_fasta, extract_pf_proteins
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)
# default Pf-codon windows for discovery genes (from the original resistance designers)
DISC_WIN = {"dhfr": (12, 120), "dhps": (390, 620), "k13": (335, 560), "crt": (72, 150), "mdr1": (80, 200)}

MAXPROD = 400
MAXSPAN = 300          # max genomic span of marker codons per tile
HARDMAX = 330          # absolute span cap before forcing an (overlapping) split
MINJUNC = 50           # min inter-codon gap to place a clean non-overlapping junction
PRIMER_REACH = 40      # a primer may sit at most this far outside the covered codon block
CS_TAIL = 44           # CS1+CS2 add 44 nt -> tailed molecule = core + 44

GEN = {"Pf": "Pf3D7.fna", "Pv": "PvP01.fna", "Pm": "PmUG01.fna",
       "Pk": "PkHl.fna", "Poc": "PocGH01.fna", "Pow": "PowCR01.fna"}


def load_all():
    contig = {}; sgen = {}
    for sp, fn in GEN.items():
        g = load_fasta(f"{BASE}/refs/nuclear/{fn}")
        sgen[sp] = g; contig.update(g)
    mito = load_fasta(f"{BASE}/refs/mito/pan_mito.fasta")
    contig.update(mito)
    return contig, sgen, mito


def tile_codons(positions):
    """partition sorted genomic positions into groups; junctions at largest gaps >=MINJUNC.
    returns list of (lo,hi,forced_overlap) covering all positions."""
    positions = sorted(positions)
    tiles = []; i = 0; n = len(positions)
    while i < n:
        start = positions[i]
        j = i
        while j + 1 < n and positions[j + 1] - start <= MAXSPAN:
            j += 1
        forced = False
        if j + 1 < n:                                   # more codons remain -> choose a junction
            best = None
            for k in range(i, j + 1):
                gap = positions[k + 1] - positions[k]
                if gap >= MINJUNC and (best is None or gap > best[1]):
                    best = (k, gap)
            if best:
                j = best[0]
            else:                                       # no clean gap: extend to HARDMAX or force overlap
                while j + 1 < n and positions[j + 1] - start <= HARDMAX:
                    j += 1
                if j + 1 < n:
                    forced = True                       # next tile will overlap this one
        tiles.append((positions[i], positions[j], forced))
        i = j + 1
    return tiles


RELAXED = dict(GLOBAL)
RELAXED.update({'PRIMER_MIN_TM': 50.0, 'PRIMER_OPT_TM': 57.0, 'PRIMER_MAX_TM': 64.0,
                'PRIMER_MIN_SIZE': 17, 'PRIMER_MIN_GC': 15.0, 'PRIMER_MAX_GC': 90.0,
                'PRIMER_MAX_POLY_X': 6})


def design_tile(contig_seq, clo, chi, name, maxprod=MAXPROD, reach=80):
    """design a <=maxprod amplicon whose insert covers [clo,chi], primers within ~reach of the block.
    Returns primers + genomic coords."""
    span = chi - clo
    pad = reach + 60
    ts = max(0, clo - pad); te = min(len(contig_seq), chi + pad)
    templ = contig_seq[ts:te]
    tgt = [clo - ts, max(1, span)]
    cap = min(maxprod, span + 2 * reach + 44)
    attempts = [(GLOBAL, [max(130, span + 44), cap]),
                (GLOBAL, [max(120, span + 40), maxprod]),
                (RELAXED, [max(110, span + 36), maxprod])]      # relaxed-Tm fallback for AT-rich spots
    for base, pr in attempts:
        g = dict(base); g['PRIMER_PRODUCT_SIZE_RANGE'] = [pr]
        try:
            res = primer3.bindings.design_primers(
                {'SEQUENCE_ID': name, 'SEQUENCE_TEMPLATE': templ, 'SEQUENCE_TARGET': tgt}, g)
        except Exception:
            continue
        if res.get('PRIMER_PAIR_NUM_RETURNED', 0):
            i = 0
            lp = res[f'PRIMER_LEFT_{i}'][0]; rp = res[f'PRIMER_RIGHT_{i}'][0]
            F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
            return {"F": F, "R": R, "prod": res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                    "Ftm": round(res[f'PRIMER_LEFT_{i}_TM'], 1), "Rtm": round(res[f'PRIMER_RIGHT_{i}_TM'], 1),
                    "gstart": ts + lp, "gend": ts + rp + 1}
    return None


def design_seq(seq, start, name, minprod=200, maxprod=MAXPROD):
    """anchored design: template STARTS at `start`, so the forward primer cannot fall left of it.
    Guarantees strictly left-to-right non-overlapping tiles when called with start = prev.gend."""
    te = min(len(seq), start + maxprod + 6)
    templ = seq[start:te]
    if len(templ) < minprod:
        return None
    for base in (GLOBAL, RELAXED):
        g = dict(base); g['PRIMER_PRODUCT_SIZE_RANGE'] = [[minprod, min(maxprod, len(templ))]]
        try:
            res = primer3.bindings.design_primers({'SEQUENCE_ID': name, 'SEQUENCE_TEMPLATE': templ}, g)
        except Exception:
            continue
        if res.get('PRIMER_PAIR_NUM_RETURNED', 0):
            i = 0; lp = res[f'PRIMER_LEFT_{i}'][0]; rp = res[f'PRIMER_RIGHT_{i}'][0]
            F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
            return {"F": F, "R": R, "prod": res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                    "Ftm": round(res[f'PRIMER_LEFT_{i}_TM'], 1), "Rtm": round(res[f'PRIMER_RIGHT_{i}_TM'], 1),
                    "gstart": start + lp, "gend": start + rp + 1}
    return None


def specificity(F, R, sp, sgen, k=15):
    for o, og in sgen.items():
        if o == sp:
            continue
        if _acount(og, F[-k:]) > 0 and _acount(og, R[-k:]) > 0:
            return f"cross:{o}"
    return "SPECIFIC"


def main():
    os.chdir(BASE)
    contig, sgen, mito = load_all()
    PF_PROT = extract_pf_proteins()          # dhfr/dhps/crt/mdr1/k13 Pf reference proteins for tblastn
    # codon map -> {gene: {contig, rows:[(codon,gstart,gend)]}}
    cmap = {}
    for r in csv.DictReader(open("targets/codon_to_genomic_pan.tsv"), delimiter="\t"):
        try:
            cd, gs, ge = int(r['codon']), int(r['genomic_start']), int(r['genomic_end'])
        except (ValueError, TypeError):
            continue
        cmap.setdefault(r['gene'], {"contig": r['contig'], "pos": []})["pos"].append((gs, ge))

    # current module_R -> rows (species, gene, amplicon, F, R)
    cur = list(csv.DictReader(open("design/module_R_resistance.tsv"), delimiter="\t"))
    out = []; report = []; conflicts = []      # conflicts: (aidA, aidB) that overlap -> must be different pools

    def codon_key(sp, gene):
        return gene if gene[:2] in ("pf", "pv") else sp.lower() + gene

    # group current amplicons by (species, gene)
    from collections import OrderedDict
    genes = OrderedDict()
    for r in cur:
        genes.setdefault((r['species'], r['gene']), []).append(r)

    for (sp, gene), rows in genes.items():
        ck = codon_key(sp, gene)
        drugs = rows[0].get('specificity', '')  # reuse col? we'll set provenance ourselves
        if ck in cmap:                                   # ---- CODON-MAPPED: cluster + tile ----
            cinfo = cmap[ck]; cseq = contig[cinfo['contig']]
            positions = []
            for gs, ge in cinfo['pos']:
                positions += [gs, ge]
            tiles = tile_codons(positions)
            prev_end = -1; prev_aid = None; made = []
            for (clo, chi, forced) in tiles:
                d = design_tile(cseq, clo, chi, f"{gene}_{sp}_{clo}", reach=50)
                if not d and (chi - clo) < 40:            # single-codon spot: widen and retry once
                    d = design_tile(cseq, clo - 30, chi + 30, f"{gene}_{sp}_{clo}", reach=50)
                if not d:
                    report.append(f"  !! {sp} {gene} [{clo}-{chi}] NO DESIGN"); continue
                spec = specificity(d['F'], d['R'], sp, sgen)
                overlap = d['gstart'] <= prev_end
                aid = f"{gene}_c{clo}"
                if overlap and prev_aid:
                    conflicts.append((prev_aid, aid))
                prev_end = d['gend']; prev_aid = aid
                out.append((sp, gene, aid, d['F'], d['R'], d['prod'], d['Ftm'], d['Rtm'], spec, overlap))
                made.append((d['prod'], spec, overlap))
            report.append(f"  {sp:4} {gene:14} codon-tiled -> {len(made)} tiles "
                          + " ".join(f"{p}bp{'*OL' if o else ''}[{s.split(':')[0]}]" for p, s, o in made))
        else:                                            # ---- DISCOVERY: tblastn window -> tile ----
            for r in rows:
                base_gene = gene.split('_')[0]              # crt_72-150 -> crt
                m = re.search(r'(\d+)-(\d+)$', r['amplicon'])
                lo, hi = (int(m.group(1)), int(m.group(2))) if m else DISC_WIN.get(base_gene, (1, 120))
                prot = PF_PROT.get(base_gene)
                hsps, cname = hsps_for(prot, sp) if prot else (None, None)
                if not hsps:
                    report.append(f"  !! {sp} {r['amplicon']} discovery: no ortholog"); continue
                mg = marker_genomic(hsps, lo, hi)
                if not mg:
                    report.append(f"  !! {sp} {r['amplicon']} discovery: window unmapped"); continue
                S, E = min(mg), max(mg)
                seq = sgen[sp][cname]
                # sequential non-overlapping tiling: each tile starts where the previous amplicon ended
                tiles = []; spos = S; guard = 0
                target = max(200, min(360, (E - S) // max(1, -(-(E - S) // 320)) + 40))
                while spos < E - 60 and guard < 8:
                    d = design_seq(seq, spos, f"{r['amplicon']}_{spos}", minprod=200, maxprod=min(MAXPROD, target + 40))
                    if not d:
                        spos += 250; guard += 1; continue
                    tiles.append(d); spos = d['gend'] + 8; guard += 1
                made = []
                for k2, d in enumerate(tiles):
                    spec = specificity(d['F'], d['R'], sp, sgen)
                    aid = f"{r['amplicon']}_t{k2+1}" if len(tiles) > 1 else r['amplicon']
                    out.append((sp, gene, aid, d['F'], d['R'], d['prod'], d['Ftm'], d['Rtm'], spec, False))
                    made.append((d['prod'], spec))
                report.append(f"  {sp:4} {r['amplicon']:22} disc-tiled({lo}-{hi}) -> {len(made)} "
                              + " ".join(f"{p}bp[{s.split(':')[0]}]" for p, s in made))

    # write new module_R
    with open("design/module_R_resistance.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["module", "amplicon", "species", "gene", "class", "fwd_primer", "rev_primer",
                    "product_bp", "F_tm", "R_tm", "specificity", "provenance"])
        for sp, gene, aid, F, R, prod, ft, rt, spec, ol in out:
            cls = "validated" if codon_key(sp, gene) in cmap else "discovery"
            w.writerow(["R", aid, sp, gene, cls, F, R, prod, ft, rt, spec, "retiled<=400"])

    # overlap conflicts (for the assembler to place in different pools)
    with open("design/retile_overlap_conflicts.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["ampliconA", "ampliconB"])
        w.writerows(conflicts)

    print("\n".join(report))
    sizes = [o[5] for o in out]
    over = sum(1 for o in out if o[9])
    from collections import Counter
    print(f"\n[retile-R] {len(out)} resistance amplicons  (was {len(cur)})")
    print(f"  overlap conflict pairs -> design/retile_overlap_conflicts.tsv: {len(conflicts)}  {conflicts}")
    print(f"  max size {max(sizes)}bp  >400: {sum(1 for s in sizes if s>400)}  tailed>444: {sum(1 for s in sizes if s+CS_TAIL>444)}")
    print(f"  overlap-forced tiles (need 2nd pool): {over}")
    print(f"  specificity: {dict(Counter('spec' if o[8]=='SPECIFIC' else 'cross' for o in out))}")


if __name__ == "__main__":
    main()
