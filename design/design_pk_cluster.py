#!/usr/bin/env python3
"""
P. knowlesi cluster-typing add-on — resolves the three subpopulations
  cluster 1 = Mf (Macaca fascicularis, Malaysian Borneo)
  cluster 2 = Mn (Macaca nemestrina, Borneo)
  cluster 3 = Peninsular Malaysia (mainland; includes the reference clone H)

Built from PRIMARY-SOURCE published markers, rebuilt as CS1/CS2-tailed SEQUENCING amplicons so they
co-run with the rest of the pan-panel. Everything is anchored by PRIMER/GENE SEQUENCE, not by the
published chromosome coordinates: those are in the PKNH (GeneDB) assembly, which is offset from the
RefSeq strain-H assembly (GCF_000006355.2) used throughout this project (e.g. NC03_2 sits ~45 kb from
its printed coordinate). Sequence anchoring makes the design assembly-version-proof.

TIER 1 — within/between-cluster diversity: the microsatellites of Divis et al. 2015 (PLoS Pathog,
  PMID 25993697, S1 Table). We take each locus' NEST-2 primer pair verbatim and keep those that still
  prime the clone-H reference (6 of 10 bind at <=3 mismatches with the correct product size; 4 have
  diverged in clone H and are dropped honestly). Trinucleotide STRs; genotyped across isolates they
  define Mf vs Mn and give relatedness/COI. NC03_2 and CD13_61 are the top differentiators (FST>0.3);
  CD08_61 sits inside the chr8 cap380 introgression window.

TIER 2 — direct cluster-diagnostic SNP amplicons (sequencing-native; we READ the fixed allele rather
  than rely on allele-specific band dropout):
    * chr1 C1A/C2F locus  — Divis 2020 fixed-difference SNP (Mf vs Mn), EID 10.3201/eid2607.190924;
      anchored on the shared reverse primer CX01AR (binds clone H 0-mismatch).
    * chr6 C2J locus      — Divis 2020 fixed-difference SNP (Mn), same paper; anchored on C206DR.
    * chr8 cap380 (PKNH_0820800) — the Mn->Mf introgression's strongest gene (~101 SNPs FST>0.4;
      Assefa 2015 PMID 26438871 / Diez Benavente 2017 10.1371/journal.pgen.1007008); anchored on the
      RefSeq gene annotation.
  Each amplicon puts the diagnostic SNP inside the read, so the sequencer calls Mf / Mn / cluster-3.

The Peninsular cluster-3 axis is additionally read from the organellar (mito cox3/cox1) ID amplicons,
whose haplotypes separate the Borneo vs Peninsular lineages (Diez Benavente 2017).
"""
import os
import sys, re, os, csv
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import primer3, in_silico_pcr as ip
from design_resistance_amplicon import GLOBAL, GENOME
from degeneracy_utils import revcomp
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)
PKG = "refs/nuclear/PkHl.fna"
CS1 = "ACACTGACGACATGGTTCTACA"; CS2 = "TACGGTAGCAGAGACTTGGTCT"


def contig_for(chrn):
    return f"NC_0119{chrn + 1:02d}.2"


TRI = re.compile(r"([ACGT]{2,4}?)\1{3,}")
# NB P. knowlesi is intrinsically microsatellite-poor (0 perfect trinucleotide runs >=6 copies
# genome-wide, vs hundreds in P. falciparum), so de-novo STR top-up is not productive; the published
# Divis loci that still prime clone-H plus the diagnostic SNP amplicons are the panel.

# TIER 1: Divis-2015 nest-2 primer pairs (F=forward-named, R=reverse-named), verbatim from S1 Table.
DIVIS = [
    ("NC03_2",   3, "GCGGGGAGGACGATAAACCATA",     "CGTCAAATGAAGAGAGCATTGCTC",    "FST>0.3 top differentiator"),
    ("CD05_06",  5, "GCCCATTGCAGCTATGCAC",        "GTTTTCGCTCCATGTTCAGCC",       "21 alleles, most polymorphic"),
    ("CD08_61",  8, "GATCAGTGGACTGGTATACACAGATA", "CCATGTAGGATGTATATTTCTTCG",    "chr8 cap380 introgression window"),
    ("NC09_1",   9, "TTGAAGCGGAATAGGGAAGGAT",     "CACACAGGTACGTGCATACATATAAG",  ""),
    ("NC10_1",  10, "CATGGCTGGTATCCCCCTGTTC",     "GGAAGGAGAGGACTAGTGCTGAAAGAG", ""),
    ("CD11_157",11, "TAGGGATCGTCAGAGGGAGG",       "GTAGCAAATATGCTCGTTCAGG",      ""),
    ("NC12_2",  12, "GTGGAAAGGAGGCAAACACACAG",    "ATGGATTCACTCACCTTGTGGTG",     ""),
    ("NC12_4",  12, "GGAGTTCTTCCGTTCCGAATG",      "AGGGGTCGTATCGTCTCGTACG",      ""),
    ("CD13_61", 13, "CCAACAACCATTCAGTTGACAAAG",   "GTTCCAATTGTTGGCCCCATTG",      "FST>0.3 top differentiator"),
    ("CD13_107",13, "GATGACCACGTTATGGATAATGTTG",  "CGTAAAACTCTGCACCTCCTTTGC",    ""),
]

# TIER 2: diagnostic loci anchored by conserved primer / gene sequence (NOT by published coordinate).
#   name, chr, forward-anchor (conserved 5' of the allele-specific fwd), reverse primer, axis, note
DIAG_PRIMER = [
    ("diag_chr1_C1A", 1, "GTTTGGTACGTTCAAGTGTTGC", "CGTCTTCCGCTTGTCGTTTTCCATGTAC",
     "Mf/Mn", "Divis2020 C1A/C2F fixed-difference locus"),
    ("diag_chr6_C2J", 6, "GATGATAATTATCTTAAAGAGCCGGA", "CAAGACATTATGAACATTGGACCGACTG",
     "Mn", "Divis2020 C2J fixed-difference locus"),
]
# cap380 anchored on RefSeq gene annotation (coordinates correct in GCF).
CAP380 = ("diag_chr8_cap380", 8, 939600, 940280, 939940, "introgression",
          "cap380 PKNH_0820800 (~101 SNPs FST>0.4)")


def locate(primer, seq, max_mm=3, seedlen=15):
    """3'-seed exact match then verify full primer <=max_mm. Returns [(start,mm)] on the + strand."""
    seed = primer[-seedlen:]; off = len(primer) - seedlen
    out = []; i = seq.find(seed)
    while i >= 0:
        ps = i - off
        if ps >= 0 and ps + len(primer) <= len(seq):
            mm = sum(1 for k in range(len(primer)) if seq[ps + k] != primer[k])
            if mm <= max_mm:
                out.append((ps, mm))
        i = seq.find(seed, i + 1)
    return out


def find_product(F, R, contig_seq, max_prod=600):
    """locate F + revcomp(R) on either strand -> smallest product (fp, rp_end, strand, insert)."""
    rr = revcomp(R)
    for strand, s in (("+", contig_seq), ("-", revcomp(contig_seq))):
        fs = locate(F, s); rs = locate(rr, s)
        best = None
        for fp, fmm in fs:
            for rp, rmm in rs:
                end = rp + len(rr); size = end - fp
                if len(F) < size <= max_prod and rp >= fp:
                    if best is None or size < best[2]:
                        best = (fp, end, size, strand, s[fp + len(F):rp], fmm, rmm)
        if best:
            return best
    return None


def _acount(g, sub):
    rc = revcomp(sub)
    return sum(c.count(sub) + c.count(rc) for c in g.values())


def uniq_and_specific(F, R, pk, others, k=15):
    fa, ra = F[-k:], R[-k:]
    nprod = 1 if (_acount(pk, fa) == 1 and _acount(pk, ra) == 1) else max(_acount(pk, fa), _acount(pk, ra))
    cross = ""
    for sp, og in others.items():
        if _acount(og, fa) > 0 and _acount(og, ra) > 0:
            cross = sp; break
    return nprod, cross


def main():
    os.chdir(BASE)
    from locate_orthologs import load_fasta
    pk = load_fasta(PKG)
    others = {}
    for sp, path in list(GENOME.items()) + [("Pf", "refs/nuclear/Pf3D7.fna")]:
        if sp != "Pk":
            try:
                others[sp] = load_fasta(path)
            except Exception:
                pass
    rows = []

    print("=== TIER 1: Divis-2015 microsatellites (primer-anchored on clone-H) ===")
    for name, chrn, F, R, note in DIVIS:
        seq = pk[contig_for(chrn)]
        prod = find_product(F, R, seq)
        if not prod:
            print(f"  {name:9} chr{chrn}: DROPPED (diverged in clone-H reference, primers do not prime)")
            continue
        fp, end, size, strand, insert, fmm, rmm = prod
        m = TRI.search(insert)
        motif = f"({m.group(1)})n" if m else "STR(short/impure in clone-H)"
        nprod, cross = uniq_and_specific(F, R, pk, others)
        spec = "Pk-specific" if not cross else f"cross:{cross}"
        rows.append((f"pkclu_{name}", "T1_microsat", chrn, motif, size, F, R, spec,
                     f"{note}; primer-anchored, {fmm}/{rmm} mm, uniq={nprod}"))
        print(f"  {name:9} chr{chrn}: {size}bp {motif} mm{fmm}/{rmm} uniq={nprod} [{spec}] {note}")

    print("=== TIER 2: cluster-diagnostic SNP amplicons ===")
    # primer-anchored diagnostics (chr1, chr6)
    for name, chrn, fanchor, rev, axis, note in DIAG_PRIMER:
        seq = pk[contig_for(chrn)]
        prod = find_product(fanchor, rev, seq, max_prod=700)
        if not prod:
            print(f"  {name}: anchor primers do not locate a product — SKIP"); continue
        fp, end, size, strand, insert, fmm, rmm = prod
        s = seq if strand == "+" else revcomp(seq)
        snp = fp + len(fanchor)                                # fixed SNP sits at the fwd-anchor 3' end
        ws = max(0, snp - 200); we = min(len(s), snp + 200)
        templ = s[ws:we]; tgt = [max(1, snp - ws - 15), 30]   # target the fixed-SNP window, not whole locus
        g = dict(GLOBAL); g['PRIMER_PRODUCT_SIZE_RANGE'] = [[220, 400]]  # platform-agnostic <=400
        res = primer3.bindings.design_primers(
            {'SEQUENCE_ID': name, 'SEQUENCE_TEMPLATE': templ, 'SEQUENCE_TARGET': tgt}, g)
        placed = None
        for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
            F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
            n, cross = uniq_and_specific(F, R, pk, others)
            if n == 1:
                placed = (F, R, res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                          round(res[f'PRIMER_LEFT_{i}_TM'], 1), round(res[f'PRIMER_RIGHT_{i}_TM'], 1),
                          "Pk-specific" if not cross else f"cross:{cross}")
                break
        if not placed:
            print(f"  {name}: no unique design — SKIP"); continue
        F, R, prod_bp, ft, rt, spec = placed
        rows.append((f"pkclu_{name}", "T2_diagnostic", chrn, axis, prod_bp, F, R, spec,
                     f"{note}; spans fixed SNP; Tm {ft}/{rt}"))
        print(f"  {name}: {prod_bp}bp axis={axis} Tm {ft}/{rt} [{spec}]  {note}")
    # gene-anchored diagnostic (chr8 cap380)
    name, chrn, ws0, we0, center, axis, note = CAP380
    seq = pk[contig_for(chrn)]
    templ = seq[ws0 - 1:we0]; tgt = [center - ws0, 30]
    g = dict(GLOBAL); g['PRIMER_PRODUCT_SIZE_RANGE'] = [[240, 400]]
    res = primer3.bindings.design_primers(
        {'SEQUENCE_ID': name, 'SEQUENCE_TEMPLATE': templ, 'SEQUENCE_TARGET': tgt}, g)
    for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
        F = res[f'PRIMER_LEFT_{i}_SEQUENCE'].upper(); R = res[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
        n, cross = uniq_and_specific(F, R, pk, others)
        if n == 1:
            rows.append((f"pkclu_{name}", "T2_diagnostic", chrn, axis,
                         res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'], F, R,
                         "Pk-specific" if not cross else f"cross:{cross}",
                         f"{note}; Tm {round(res[f'PRIMER_LEFT_{i}_TM'],1)}/{round(res[f'PRIMER_RIGHT_{i}_TM'],1)}"))
            print(f"  {name}: {res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}bp axis={axis} [{'Pk-specific' if not cross else 'cross:'+cross}]  {note}")
            break

    out = f"{BASE}/design/module_V_pk_cluster.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["amplicon", "tier", "chr", "motif_or_axis", "product_bp",
                    "fwd_primer", "rev_primer", "specificity", "evidence"])
        w.writerows(rows)
    from collections import Counter
    print(f"\n[pk-cluster] {len(rows)} amplicons -> {out}")
    print("  tiers:", dict(Counter(r[1] for r in rows)),
          "| Pk-specific:", sum(1 for r in rows if r[7] == "Pk-specific"))


if __name__ == "__main__":
    main()
