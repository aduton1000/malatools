#!/usr/bin/env python3
"""
build_manifest.py — one-time index step.

Locate every panel amplicon on the combined pan-reference (refs/pan_reference.fa.gz) and write
analysis/manifest/amplicons.tsv + amplicons.bed, used by every downstream analysis stage to pile up reads
at the right coordinates. Robust to contig-naming: we extract each amplicon's expected sequence from the
species genome, then minimap2 it to the pan-reference to get exact coordinates on the pan contig.

Columns: amplicon  module  species  gene  role  pan_contig  start  end  strand  product_bp
For the pan-species cox3 species-ID amplicon, one row per species (amplicon = panPlas_cox3_S1__<sp>).
"""
import sys, os, csv, subprocess, gzip
HERE = os.path.dirname(os.path.abspath(__file__))
DESIGN = os.path.abspath(os.path.join(HERE, "../../design"))
BASE = os.path.abspath(os.path.join(HERE, "../.."))
sys.path.insert(0, DESIGN)
import in_silico_pcr as ip
from degeneracy_utils import revcomp
from design_pk_cluster import find_product
OUT = os.path.join(BASE, "analysis/manifest")
os.makedirs(OUT, exist_ok=True)

GEN = {"Pf": "Pf3D7.fna", "Pv": "PvP01.fna", "Pm": "PmUG01.fna",
       "Pk": "PkHl.fna", "Poc": "PocGH01.fna", "Pow": "PowCR01.fna"}
MITO_ACC = {"Pf": "NC_002375.1", "Pv": "NC_007243.1", "Poc": "CM190334.1",
            "Pow": "CM190336.1", "Pm": "LT594637.1", "Pk": "NC_007232.1"}
ALL_SP = list(GEN)


def load_genomes():
    G = {}
    for sp, fn in GEN.items():
        G[sp] = ip.load_fasta(f"{BASE}/refs/nuclear/{fn}")
    mito = ip.load_fasta(f"{BASE}/refs/mito/pan_mito.fasta")
    for sp, acc in MITO_ACC.items():
        G[sp]["__mito__:" + acc] = mito[acc]      # attach species mito under a tagged key
    return G


def mito_key(sp):
    return "__mito__:" + MITO_ACC[sp]


def contigs_for(sp, gene, role, G):
    """restrict the search space: cox3/cytb -> mito only; else nuclear (avoids scanning 654/780
    fragmented Poc/Pow contigs with the slow IUPAC fallback)."""
    g = (gene or "").lower(); mk = mito_key(sp)
    if role == "species_id" or "cox3" in g or "cox1" in g or "cytb" in g:
        return {mk: G[sp][mk]}
    return {k: v for k, v in G[sp].items() if not k.startswith("__mito__")}


def extract(F, R, contigs, allow_iupac):
    """locate F/R in {name:seq}; fast exact-seed first, IUPAC/mismatch fallback only when allow_iupac
    (small contig sets, e.g. mito) to avoid O(contigs) pure-python isPCR on fragmented assemblies."""
    for name, seq in contigs.items():
        fp = find_product(F, R, seq, max_prod=1600)
        if fp:
            strand = fp[3]; s = seq if strand == "+" else revcomp(seq)
            return s[fp[0]:fp[1]]
    if allow_iupac:
        for name, seq in contigs.items():
            hits = ip.ispcr_bothstrands(F, R, seq, max_mm=4, clean3=2, min_prod=120, max_prod=1600)
            if hits:
                strand, p = min(hits, key=lambda x: x[1]["size"])
                s = seq if strand == "+" else revcomp(seq)
                return s[p["start"]:p["end"]]
    return None


def collect_amplicons():
    """yield (amplicon_id, module, species, gene, role, F, R) from all module tables."""
    def rows(f):
        p = f"{DESIGN}/{f}"
        return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []
    for r in rows("module_S_species_id.tsv"):
        if r["fwd_primer"] and r["rev_primer"] and not r.get("role", "").startswith("OPTIONAL"):
            for sp in ALL_SP:                          # cox3 amplifies all 6 taxa
                yield (f"{r['amplicon']}__{sp}", "S", sp, "cox3", "species_id", r["fwd_primer"], r["rev_primer"])
    for r in rows("module_R_resistance.tsv"):
        if r.get("fwd_primer") and r.get("rev_primer"):
            yield (r["amplicon"], "R", r["species"], r["gene"], "resistance", r["fwd_primer"], r["rev_primer"])
    for r in rows("module_D_diagnostic_cnv.tsv"):
        if r.get("fwd_primer") and r.get("rev_primer"):
            yield (r["amplicon"], "D", r["species"], r["gene"], r["role"], r["fwd_primer"], r["rev_primer"])
    for r in rows("module_V_amplicons.tsv"):
        if r.get("fwd_primer") and r.get("rev_primer"):
            yield (r["amplicon"], "V", r["species"], "microhap", "diversity", r["fwd_primer"], r["rev_primer"])
    for r in rows("module_V_newsp_diversity.tsv"):
        if r.get("fwd_primer") and r.get("rev_primer"):
            yield (r["amplicon"], "V", r["species"], r["gene"], "diversity", r["fwd_primer"], r["rev_primer"])
    for r in rows("module_V_pk_cluster.tsv"):
        if r.get("fwd_primer") and r.get("rev_primer"):
            role = "pk_cluster_" + ("diag" if r["tier"].startswith("T2") else "microsat")
            yield (r["amplicon"], "V", "Pk", r.get("motif_or_axis", "cluster"), role, r["fwd_primer"], r["rev_primer"])


def main():
    G = load_genomes()
    recs = list(collect_amplicons())
    print(f"[manifest] {len(recs)} amplicon records; extracting sequences...")
    fa = os.path.join(OUT, "amplicons.fa")
    meta = {}
    n_ext = 0
    with open(fa, "w") as fh:
        for aid, mod, sp, gene, role, F, R in recs:
            contigs = contigs_for(sp, gene, role, G)
            prod = extract(F, R, contigs, allow_iupac=(len(contigs) <= 3))
            if not prod:
                print(f"  !! {aid} ({sp}): not extracted"); continue
            fh.write(f">{aid}\n{prod}\n")
            meta[aid] = (mod, sp, gene, role, len(prod))
            n_ext += 1
    print(f"[manifest] extracted {n_ext}/{len(recs)} amplicon sequences -> {fa}")

    # map to pan-reference for coordinates
    ref = f"{BASE}/refs/pan_reference.fa.gz"
    paf = os.path.join(OUT, "amplicons.paf")
    with open(paf, "w") as fh:
        subprocess.run(["minimap2", "-cx", "sr", "--secondary=no", "-t", "4", ref, fa],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)
    best = {}
    for line in open(paf):
        f = line.split("\t")
        q, ql, qs, qe, strand, t, tl, ts, te, nm, bl = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
        score = nm
        if q not in best or score > best[q][0]:
            best[q] = (score, t, ts, te, strand)
    mtsv = os.path.join(OUT, "amplicons.tsv"); bed = os.path.join(OUT, "amplicons.bed")
    nmap = 0
    with open(mtsv, "w", newline="") as fh, open(bed, "w") as bh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["amplicon", "module", "species", "gene", "role", "pan_contig", "start", "end", "strand", "product_bp"])
        for aid, (mod, sp, gene, role, plen) in meta.items():
            if aid not in best:
                print(f"  !! {aid}: no pan-reference hit"); continue
            sc, t, ts, te, strand = best[aid]
            w.writerow([aid, mod, sp, gene, role, t, ts, te, strand, plen])
            bh.write(f"{t}\t{ts}\t{te}\t{aid}\t0\t{strand}\n")
            nmap += 1
    print(f"[manifest] mapped {nmap} amplicons -> {mtsv} + {bed}")
    from collections import Counter
    print("  by module:", dict(Counter(v[0] for v in meta.values())))


if __name__ == "__main__":
    main()
