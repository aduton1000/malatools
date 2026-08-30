#!/usr/bin/env python3
"""
validate.py — run the pipeline output against the simulated ground truth and report PASS/FAIL for every
assertion across all modules (species ID incl. mixed, resistance incl. reference-is-mutant polarity, CNV,
gene deletion, COI, Pk cluster). Exit non-zero if any assertion fails.

Usage: validate.py <sim_dir> <out_dir>
"""
import sys, os, csv
from collections import defaultdict


def read_tsv(path, header=True):
    if not os.path.exists(path):
        return []
    rows = list(csv.reader(open(path), delimiter="\t"))
    return rows[1:] if header and rows else rows


def main():
    sim, out = sys.argv[1], sys.argv[2]
    T = {k: read_tsv(f"{sim}/truth/{k}.tsv", header=False) for k in
         ["species", "resistance", "cnv", "deletion", "coi", "cluster", "related"]}
    samples = sorted({r[0] for r in T["species"]})
    P = []; F = []

    def check(cond, msg):
        (P if cond else F).append(msg)

    for s in samples:
        od = f"{out}/{s}"
        # ---- SPECIES ----
        truth_present = {r[1] for r in T["species"] if r[0] == s and r[2] == "present"}
        got = read_tsv(f"{od}/{s}.species.tsv")
        called = {r[1] for r in got if len(r) > 4 and r[4] == "present"}
        check(called == truth_present,
              f"[{s}] species: expected {sorted(truth_present) or 'none'} got {sorted(called) or 'none'}")

        # ---- RESISTANCE (incl. polarity) ----
        res = {(r[1], r[2], r[3]): r for r in read_tsv(f"{od}/{s}.resistance.tsv")}   # (sp,gene,codon)->row
        for r in [x for x in T["resistance"] if x[0] == s]:
            _, sp, gene, codon, taa = r
            row = res.get((sp, gene, codon))
            if not row:
                F.append(f"[{s}] resistance {sp} {gene} {codon}: no output row"); continue
            obs, wt, call = row[6], row[5], row[9]
            if wt == "?":
                check(obs == taa, f"[{s}] {sp} {gene}{codon}: genotype expected {taa} got {obs}")
            elif taa == wt:      # sensitive truth (may be reference-is-mutant locus)
                check(call == "sensitive", f"[{s}] {sp} {gene}{codon}: expected SENSITIVE (truth {taa}=wt) got {call}/{obs}")
            else:                # resistant truth
                check(call == "RESISTANT" and obs == taa,
                      f"[{s}] {sp} {gene}{codon}: expected RESISTANT {taa} got {call}/{obs}")

        # ---- CNV ----
        cnv = {r[2]: r for r in read_tsv(f"{od}/{s}.cnv.tsv")}
        for r in [x for x in T["cnv"] if x[0] == s]:
            _, sp, gene, cn = r
            row = cnv.get(gene)
            check(row is not None and "AMPLIFIED" in row[7],
                  f"[{s}] CNV {gene}: expected AMPLIFIED (truth {cn} copies) got {row[7] if row else 'missing'}")

        # ---- DELETION ----
        for r in [x for x in T["deletion"] if x[0] == s]:
            _, sp, gene, exp = r
            row = cnv.get(gene)
            check(row is not None and (("DELETED" in row[7]) == (exp == "DELETED")),
                  f"[{s}] deletion {gene}: expected {exp} got {row[7] if row else 'missing'}")

        # ---- COI ----
        div = {r[1]: r for r in read_tsv(f"{od}/{s}.diversity.tsv")}
        for r in [x for x in T["coi"] if x[0] == s]:
            _, sp, k = r; k = int(k)
            row = div.get(sp)
            if not row:
                if k > 1:
                    F.append(f"[{s}] COI {sp}: expected ~{k} but no diversity output")
                continue
            coi = int(row[6])
            check(abs(coi - k) <= 1, f"[{s}] COI {sp}: expected ~{k} got {coi}")

        # ---- PK CLUSTER ----
        for r in [x for x in T["cluster"] if x[0] == s]:
            exp = r[1]
            rows = read_tsv(f"{od}/{s}.cluster.tsv")
            callrow = next((x for x in rows if x[1] == "CALL"), None)
            check(callrow is not None and callrow[5].startswith(exp),
                  f"[{s}] Pk cluster: expected {exp} got {callrow[5] if callrow else 'missing'}")

    # ---- COHORT: COI-MLE + relatedness (from relatedness.py) ----
    import glob as _glob
    coif = _glob.glob(f"{out}/*.coi.tsv")
    if coif:
        coi_mle = {(r[0], r[1]): int(r[2]) for r in read_tsv(coif[0])}
        for r in T["coi"]:
            s, sp, k = r[0], r[1], int(r[2])
            m = coi_mle.get((s, sp))
            if m is not None:
                check(abs(m - k) <= 1, f"[{s}] COI-MLE {sp}: expected ~{k} got {m}")
    relf = _glob.glob(f"{out}/*.relatedness.tsv")
    if relf and T["related"]:
        rel = {frozenset((r[0], r[1])): float(r[4]) for r in read_tsv(relf[0])}
        truthset = {frozenset(p) for p in T["related"]}
        for p in T["related"]:
            v = rel.get(frozenset(p))
            check(v is not None and v >= 0.7, f"related pair {p[0]}~{p[1]}: expected r>=0.7 got {v}")
        maxun = max([v for k, v in rel.items() if k not in truthset], default=0.0)
        check(maxun <= 0.3, f"relatedness specificity: unrelated pairs max r should be <=0.3, got {maxun}")

    print("=" * 72)
    print(f"VALIDATION: {len(P)} PASS, {len(F)} FAIL  ({len(samples)} samples)")
    print("=" * 72)
    for m in F:
        print("  FAIL ", m)
    if not F:
        print("  ✓ every module recovered ground truth on all samples")
    print()
    # concise per-category tally
    def catkey(m):
        body = m.split("]", 1)[1] if "]" in m else m
        return body.strip().split(" ")[0].split(":")[0]
    cat = defaultdict(lambda: [0, 0])
    for m in P:
        cat[catkey(m)][0] += 1
    for m in F:
        cat[catkey(m)][1] += 1
    for k, (p, f) in sorted(cat.items()):
        print(f"  {k:14} {p} pass / {f} fail")
    sys.exit(1 if F else 0)


if __name__ == "__main__":
    main()
