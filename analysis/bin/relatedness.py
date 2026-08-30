#!/usr/bin/env python3
"""
relatedness.py — cohort-level complexity-of-infection (COI) and pairwise relatedness from the per-sample
microhaplotype allele calls (<sample>.alleles.tsv from diversity.py). Native, dependency-free estimators:

  * COI (maximum likelihood): for each sample & species, the number of clones m that best explains the
    observed set of microhaplotypes at each locus given population allele frequencies — via the
    inclusion-exclusion probability that m draws cover exactly the observed allele set. This upgrades the
    per-sample heuristic COI to a proper MLE (the same idea as THE REAL McCOIL's categorical method).
  * Relatedness (pairwise r): a population-allele-frequency-adjusted identity-by-state moment estimator on
    the dominant-clone microhaplotypes — r ≈ (observed matches − chance matches) / (1 − chance matches).

Also EXPORTS the microhaplotype tables in the input formats of the canonical tools, so publication-grade
runs are one command away:  THE REAL McCOIL (categorical matrix) and Dcifer (long allele table).

Usage: relatedness.py <results_dir> <cohort_name> [--present 0.1] [--min-loci 5] [--mmax 6] [--r-report 0.25]
"""
import os, sys, csv, glob, math, argparse
from itertools import combinations
from collections import defaultdict


def load_alleles(results, present):
    """-> data[species][sample][locus] = set(present alleles); majors[species][sample][locus]=major allele"""
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    major = defaultdict(lambda: defaultdict(dict))
    for f in glob.glob(f"{results}/*/*.alleles.tsv"):
        rows = list(csv.DictReader(open(f), delimiter="\t"))
        best = {}
        for r in rows:
            s, sp, loc, al, fr = r["sample"], r["species"], r["locus"], r["allele"], float(r["read_fraction"])
            if fr >= present:
                data[sp][s][loc].add(al)
            if (sp, s, loc) not in best or fr > best[(sp, s, loc)][1]:
                best[(sp, s, loc)] = (al, fr)
        for (sp, s, loc), (al, fr) in best.items():
            major[sp][s][loc] = al
    return data, major


def load_consensus(results):
    cons = defaultdict(lambda: defaultdict(dict))       # species -> sample -> locus -> consensus seq
    for f in glob.glob(f"{results}/*/*.consensus.tsv"):
        for r in csv.DictReader(open(f), delimiter="\t"):
            cons[r["species"]][r["sample"]][r["locus"]] = r["consensus"]
    return cons


def build_microhaps(cons_sp):
    """from per-sample per-locus consensus, define each sample's microhaplotype allele = its bases at the
    positions that vary ACROSS the cohort (population-informative sites). -> majors[sample][locus]=allele."""
    loci = defaultdict(dict)
    for s, ld in cons_sp.items():
        for loc, seq in ld.items():
            loci[loc][s] = seq
    majors = defaultdict(dict)
    for loc, seqs in loci.items():
        if len(seqs) < 2:
            continue
        L = min(len(s) for s in seqs.values())
        vpos = [j for j in range(L) if len({seqs[s][j] for s in seqs if seqs[s][j] != "N"}) > 1]
        if not vpos:
            continue                                     # monomorphic across cohort -> uninformative
        for s, seq in seqs.items():
            majors[s][loc] = "".join(seq[j] for j in vpos)
    return majors


def pop_freq(sample_locus_alleles):
    """population allele frequency per locus (presence-based, normalised to sum 1)."""
    counts = defaultdict(lambda: defaultdict(int))
    for s, loci in sample_locus_alleles.items():
        for loc, alset in loci.items():
            for a in alset:
                counts[loc][a] += 1
    freq = {}
    for loc, c in counts.items():
        tot = sum(c.values())
        freq[loc] = {a: n / tot for a, n in c.items()} if tot else {}
    return freq


def major_freq(major_by_sample):
    counts = defaultdict(lambda: defaultdict(int))
    for s, loci in major_by_sample.items():
        for loc, a in loci.items():
            counts[loc][a] += 1
    freq = {}
    for loc, c in counts.items():
        tot = sum(c.values()); freq[loc] = {a: n / tot for a, n in c.items()}
    return freq


def coi_mle(sample_loci, pop, mmax, eps=0.05):
    """MLE COI. `eps` reserves probability mass for unseen alleles so that P(observed set) < 1 and the
    likelihood has an interior maximum (otherwise, with few samples / low-diversity loci, m runs to mmax)."""
    best_m, best_ll = 1, -1e18
    for m in range(1, mmax + 1):
        ll, valid = 0.0, True
        for loc, S in sample_loci.items():
            k = len(S)
            if k == 0:
                continue
            if m < k:
                valid = False; break
            p = pop.get(loc, {})
            Sl = list(S); P = 0.0
            for r in range(1, k + 1):
                for T in combinations(Sl, r):
                    pt = (1 - eps) * sum(p.get(a, 1e-6) for a in T)   # reserve eps for unseen alleles
                    P += ((-1) ** (k - r)) * (pt ** m)
            ll += math.log(max(P, 1e-12))
        if valid and ll > best_ll:
            best_ll, best_m = ll, m
    return best_m


def relatedness(gi, gj, mfreq, min_loci):
    shared = [l for l in gi if l in gj]
    if len(shared) < min_loci:
        return None, len(shared)
    Sobs = sum(1 for l in shared if gi[l] == gj[l])
    Sexp = sum(sum(f * f for f in mfreq.get(l, {}).values()) for l in shared)
    n = len(shared)
    if n - Sexp <= 0:
        return 0.0, n
    return max(0.0, min(1.0, (Sobs - Sexp) / (n - Sexp))), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results"); ap.add_argument("cohort")
    ap.add_argument("--present", type=float, default=0.1); ap.add_argument("--min-loci", type=int, default=5)
    ap.add_argument("--mmax", type=int, default=6); ap.add_argument("--r-report", type=float, default=0.25)
    a = ap.parse_args()
    data, _ = load_alleles(a.results, a.present)
    cons = load_consensus(a.results)
    if not data:
        sys.exit("[relatedness] no *.alleles.tsv found (run diversity first)")

    coi_rows = []; rel_rows = []; pairs = []
    dcifer = [["sample.id", "locus", "allele", "species"]]
    for sp in sorted(data):
        pf = pop_freq(data[sp])
        major = {sp: build_microhaps(cons[sp])}          # population-informative microhaplotypes
        mf = major_freq(major[sp])
        samples = sorted(data[sp])
        # COI
        for s in samples:
            m = coi_mle(data[sp][s], pf, a.mmax)
            nloci = len(data[sp][s])
            coi_rows.append([s, sp, m, nloci])
            for loc, alset in data[sp][s].items():
                for al in alset:
                    dcifer.append([s, loc, al, sp])
        # relatedness (all pairs)
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                r, n = relatedness(major[sp].get(samples[i], {}), major[sp].get(samples[j], {}), mf, a.min_loci)
                if r is None:
                    continue
                rel_rows.append([samples[i], samples[j], sp, n, round(r, 3)])
                if r >= a.r_report:
                    pairs.append([samples[i], samples[j], sp, n, round(r, 3)])

    os.makedirs(a.results, exist_ok=True)
    with open(f"{a.results}/{a.cohort}.coi.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["sample", "species", "COI_mle", "loci_typed"]); w.writerows(coi_rows)
    with open(f"{a.results}/{a.cohort}.relatedness.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["sample_i", "sample_j", "species", "shared_loci", "relatedness_r"]); w.writerows(rel_rows)
    with open(f"{a.results}/{a.cohort}.related_pairs.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["sample_i", "sample_j", "species", "shared_loci", "relatedness_r"]); w.writerows(sorted(pairs, key=lambda x: -x[4]))

    # ---- exports for the canonical tools ----
    exp = f"{a.results}/export"; os.makedirs(exp, exist_ok=True)
    with open(f"{exp}/dcifer_input.tsv", "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(dcifer)
    # THE REAL McCOIL categorical: samples x loci, 0=single allele, 0.5=het(>1 allele), -1=missing
    for sp in sorted(data):
        loci = sorted({l for s in data[sp] for l in data[sp][s]})
        with open(f"{exp}/realmccoil_{sp}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t"); w.writerow(["sample"] + loci)
            for s in sorted(data[sp]):
                row = [s]
                for l in loci:
                    al = data[sp][s].get(l)
                    row.append(-1 if not al else (0.5 if len(al) > 1 else 0))
                w.writerow(row)

    sys.stderr.write(f"[relatedness] {len(coi_rows)} sample COIs, {len(rel_rows)} pairs "
                     f"({len(pairs)} with r>={a.r_report}); exports in {exp}\n")
    for p in sorted(pairs, key=lambda x: -x[4])[:5]:
        sys.stderr.write(f"    related: {p[0]} ~ {p[1]} ({p[2]}) r={p[4]} over {p[3]} loci\n")


if __name__ == "__main__":
    main()
