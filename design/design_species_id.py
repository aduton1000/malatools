#!/usr/bin/env python3
"""
Design a pan-Plasmodium species-ID amplicon on a mito-gene MSA.

Strategy: pick a forward conserved window (5' third) and a reverse conserved window (3' third) with a
product in the target size range; build DEGENERATE primers (IUPAC consensus, capped fold, clean 3');
then verify the INSERT between primers discriminates all six taxa (every pairwise distance > 0), with
special attention to the hardest pair (P. ovale curtisi vs wallikeri). Score = maximise the minimum
pairwise discrimination while minimising primer degeneracy and balancing Tm.

Usage: design_species_id.py <gene.aln.fasta> <min_prod> <max_prod>
"""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import degeneracy_utils as du

ORDER = ["Pf", "Pv", "Poc", "Pow", "Pm", "Pk"]

def read_aln(p):
    d = {}
    for r in open(p).read().split(">"):
        if r.strip():
            h, s = r.split("\n", 1); d[h.strip()] = s.replace("\n", "").upper()
    return d

def windows(aln_list, lo, hi, plen_range=(20, 27), max_deg=16):
    """conserved candidate windows across plen options, indexed by start."""
    out = []
    for plen in range(plen_range[0], plen_range[1] + 1):
        for w in du.primer_windows(aln_list, plen=plen, max_deg=max_deg, tail3_clean=5):
            if lo <= w["start"] <= hi:
                w["plen"] = plen; out.append(w)
    return out

def insert_discrimination(seqs):
    """seqs: list of aligned insert strings (same length). Return (min_pairwise, worst_pair, matrix)."""
    n = len(seqs); mat = [[0] * n for _ in range(n)]
    mn = 10**9; worst = None
    for i in range(n):
        for j in range(i + 1, n):
            d = sum(1 for a, b in zip(seqs[i], seqs[j]) if a != b)
            mat[i][j] = mat[j][i] = d
            if d < mn: mn = d; worst = (ORDER[i], ORDER[j])
    return mn, worst, mat

def tm_ok(f, r, min_mean=48, max_spread=8):
    return (f["tm_mean"] >= min_mean and r["tm_mean"] >= min_mean
            and abs(f["tm_mean"] - r["tm_mean"]) <= max_spread)

def main():
    aln_path, minp, maxp = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    d = read_aln(aln_path)
    seqs = [d[k] for k in ORDER if k in d]; present = [k for k in ORDER if k in d]
    L = len(seqs[0])
    fwd = windows(seqs, 0, int(L * 0.45))
    rev = windows(seqs, int(L * 0.35), L - 20)
    best = None
    for f in fwd:
        fe = f["start"] + f["plen"]
        for r in rev:
            prod = (r["start"] + r["plen"]) - f["start"]
            if not (minp <= prod <= maxp): continue
            if r["start"] < fe + 40: continue                 # need an insert
            if not tm_ok(f, r): continue
            ins = [s[fe:r["start"]] for s in seqs]
            if any(("-" in x) for x in ins):                  # avoid gappy inserts near primers is fine; gaps counted
                pass
            mn, worst, mat = insert_discrimination(ins)
            if mn <= 0: continue                              # some taxa identical -> reject
            score = mn - 0.5 * (du.fold_degeneracy(f["primer"]) + du.fold_degeneracy(r["primer"]))
            cand = {"f": f, "r": r, "prod": prod, "min_disc": mn, "worst": worst,
                    "mat": mat, "score": score, "ins": ins}
            if best is None or score > best["score"]:
                best = cand
    if not best:
        print("No amplicon found in range"); return
    f, r = best["f"], best["r"]
    fwd_primer = f["primer"]; rev_primer = du.revcomp(r["primer"])
    print(f"# species-ID amplicon on {aln_path.split('/')[-1]}  (taxa: {','.join(present)})")
    print(f"product ~{best['prod']} bp | min pairwise discrimination {best['min_disc']} SNP "
          f"(hardest pair {best['worst'][0]} vs {best['worst'][1]})")
    print(f"FWD  {fwd_primer}  (len {f['plen']}, fold {du.fold_degeneracy(fwd_primer)}, "
          f"Tm {f['tm_min']}-{f['tm_max']}, GC {f['gc_min']}-{f['gc_max']})")
    print(f"REV  {rev_primer}  (len {r['plen']}, fold {du.fold_degeneracy(r['primer'])}, "
          f"Tm {r['tm_min']}-{r['tm_max']}, GC {r['gc_min']}-{r['gc_max']})")
    print("\npairwise insert distances (SNPs):")
    print("      " + "  ".join(f"{k:>4}" for k in present))
    for i, k in enumerate(present):
        print(f"  {k:>4} " + "  ".join(f"{best['mat'][i][j]:>4}" for j in range(len(present))))

if __name__ == "__main__":
    main()
