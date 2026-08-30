#!/usr/bin/env python3
"""
Panel assembly & QC:
  1. collect all amplicon primers (Modules S, R, D), tail with CS1/CS2
  2. HUMAN off-target: in-silico PCR every pair against human mito + human DHFR
  3. cross-dimer screen (primer3 end-stability) on the tailed primers
  4. pool assignment (greedy graph colouring) if dimers force splitting
Outputs: design/panel_tailed.tsv, design/panel_qc_dimers.tsv, design/panel_pools.tsv
"""
import os
import sys, csv, os
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import primer3, in_silico_pcr as ip
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)
os.chdir(BASE)
CS1 = "ACACTGACGACATGGTTCTACA"; CS2 = "TACGGTAGCAGAGACTTGGTCT"
DIMER_TM = 40.0     # end-stability Tm above which a cross-primer dimer is flagged

def collect():
    amps = []
    for r in csv.DictReader(open("design/module_S_species_id.tsv"), delimiter="\t"):
        if r["fwd_primer"] and r["rev_primer"] and not r.get("role", "").startswith("OPTIONAL"):
            amps.append((r["amplicon"], "S", r["fwd_primer"], r["rev_primer"]))   # OPTIONAL-ONT excluded from universal pool
    for r in csv.DictReader(open("design/module_R_resistance.tsv"), delimiter="\t"):
        if r["fwd_primer"] and r["rev_primer"]:
            amps.append((r["amplicon"], "R", r["fwd_primer"], r["rev_primer"]))
    for r in csv.DictReader(open("design/module_D_diagnostic_cnv.tsv"), delimiter="\t"):
        if r.get("fwd_primer") and r.get("rev_primer"):
            amps.append((r["amplicon"], "D", r["fwd_primer"], r["rev_primer"]))
    for r in csv.DictReader(open("design/module_V_amplicons.tsv"), delimiter="\t"):
        if r.get("fwd_primer") and r.get("rev_primer"):
            amps.append((r["amplicon"], "V", r["fwd_primer"], r["rev_primer"]))
    # new-species diversity (discovery-grade, homology-selected antigen loci)
    nsp = "design/module_V_newsp_diversity.tsv"
    if os.path.exists(nsp):
        for r in csv.DictReader(open(nsp), delimiter="\t"):
            if r.get("fwd_primer") and r.get("rev_primer"):
                amps.append((r["amplicon"], "V", r["fwd_primer"], r["rev_primer"]))
    # P. knowlesi cluster-typing add-on (Divis microsatellites + diagnostic SNP amplicons)
    pkc = "design/module_V_pk_cluster.tsv"
    if os.path.exists(pkc):
        for r in csv.DictReader(open(pkc), delimiter="\t"):
            if r.get("fwd_primer") and r.get("rev_primer"):
                amps.append((r["amplicon"], "V", r["fwd_primer"], r["rev_primer"]))
    # de-dup by amplicon id
    seen = {}; out = []
    for a in amps:
        if a[0] not in seen: seen[a[0]] = 1; out.append(a)
    return out

def main():
    amps = collect()
    # 1. tailed table
    with open("design/panel_tailed.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["amplicon", "module", "core_F", "core_R", "tailed_F(CS1+F)", "tailed_R(CS2+R)"])
        for aid, mod, F, R in amps:
            w.writerow([aid, mod, F, R, CS1 + F, CS2 + R])
    print(f"[assemble] {len(amps)} amplicons tailed -> panel_tailed.tsv")

    # 2. human off-target
    human = ip.load_fasta("refs/human_offtarget.fasta")
    hits = []
    for aid, mod, F, R in amps:
        for hn, hs in human.items():
            for strand, p in ip.ispcr_bothstrands(F, R, hs, max_mm=5, clean3=3, min_prod=50, max_prod=3000):
                hits.append((aid, hn, p["size"]))
    print(f"[assemble] human off-target products: {len(hits)}"
          + (" -> " + "; ".join(f"{a}:{h}:{s}bp" for a, h, s in hits) if hits else " (none)"))

    # 3. cross-dimer screen on tailed primers (flag 3'-involving dimers between DIFFERENT amplicons)
    prim = []
    for aid, mod, F, R in amps:
        prim.append((aid, "F", CS1 + F)); prim.append((aid, "R", CS2 + R))
    flags = []
    incompat = {aid: set() for aid, _, _, _ in amps}
    for i in range(len(prim)):
        a1, d1, s1 = prim[i]
        for j in range(i + 1, len(prim)):
            a2, d2, s2 = prim[j]
            if a1 == a2:
                continue
            tm = primer3.bindings.calc_end_stability(s1, s2).tm
            if tm > DIMER_TM:
                flags.append((round(tm, 1), f"{a1}.{d1}", f"{a2}.{d2}"))
                incompat[a1].add(a2); incompat[a2].add(a1)
    # overlap edges: tiled amplicons of the same locus must be in different pools
    import re as _re
    def base(aid):
        b = _re.sub(r"_(A|S|D)\d+$", "", aid)            # pfcrt_A1/pfcrt_A2 -> pfcrt ; cox3_S1 -> cox3
        b = _re.sub(r"_\d+-\d+$", "|W", b)               # keep mdr1 windows distinct (non-overlap)
        return b
    from collections import defaultdict as _dd
    grp = _dd(list)
    for aid, mod, F, R in amps: grp[base(aid)].append(aid)
    for b, members in grp.items():
        if len(members) > 1 and not b.endswith("|W"):
            for x in members:
                for y in members:
                    if x != y: incompat[x].add(y)
    # explicit overlap conflicts from the resistance re-tiler (genuinely-overlapping tiles, e.g. pfk13)
    ov = "design/retile_overlap_conflicts.tsv"
    if os.path.exists(ov):
        for r in csv.DictReader(open(ov), delimiter="\t"):
            a, b2 = r["ampliconA"], r["ampliconB"]
            if a in incompat and b2 in incompat:
                incompat[a].add(b2); incompat[b2].add(a)
    flags.sort(reverse=True)
    with open("design/panel_qc_dimers.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["end_stability_Tm", "primerA", "primerB"])
        w.writerows(flags)
    print(f"[assemble] cross-primer dimers > {DIMER_TM}C: {len(flags)}  (top: "
          + ("; ".join(f'{t}C {a}x{b}' for t, a, b in flags[:5]) if flags else "none") + ")")

    # 4. pool assignment (greedy colouring by incompatibility)
    order = sorted(incompat, key=lambda a: -len(incompat[a]))
    pool = {}
    for a in order:
        used = {pool[n] for n in incompat[a] if n in pool}
        k = 1
        while k in used: k += 1
        pool[a] = k
    npools = max(pool.values()) if pool else 1
    with open("design/panel_pools.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["amplicon", "pool"])
        for aid, mod, F, R in amps: w.writerow([aid, pool.get(aid, 1)])
    from collections import Counter
    print(f"[assemble] pools needed: {npools}  {dict(Counter(pool.values()))}")

if __name__ == "__main__":
    main()
