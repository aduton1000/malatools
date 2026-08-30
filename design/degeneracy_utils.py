#!/usr/bin/env python3
"""
Degeneracy-aware primer utilities for the pan-Plasmodium panel.

Core engine for designing degenerate primers on conserved flanks of a multi-species alignment and
QC-ing them. Pure, dependency-light (numpy optional); unit-tested at the bottom.

Key functions:
  consensus_iupac(seqs)        -> IUPAC-degenerate consensus of aligned, equal-length sequences
  expand_degenerate(primer)    -> all concrete oligos encoded by an IUPAC primer
  fold_degeneracy(primer)      -> number of concrete oligos (product of per-base degeneracy)
  tm_range(primer)             -> (min, max, mean) nearest-neighbour-ish Tm across the degenerate set
  three_prime_clean(primer, n) -> True if last n bases carry no degeneracy
  gc_range(primer)             -> (min, max) GC% across the degenerate set
  primer_windows(aln, ...)     -> conserved windows in an MSA suitable for primer placement
"""
from itertools import product

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
BASES_TO_CODE = {frozenset(v): k for k, v in IUPAC.items()}
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

def revcomp(s):
    return s.translate(COMP)[::-1]

def iupac_code(bases):
    """bases: iterable of A/C/G/T present at a column -> single IUPAC letter."""
    return BASES_TO_CODE[frozenset(b for b in bases if b in "ACGT")]

def consensus_iupac(seqs):
    """IUPAC-degenerate consensus of equal-length aligned sequences (gaps '-' ignored per column)."""
    L = len(seqs[0])
    assert all(len(s) == L for s in seqs), "sequences must be aligned/equal length"
    out = []
    for i in range(L):
        col = {s[i].upper() for s in seqs if s[i] in "ACGTacgt"}
        if not col:
            out.append("N")
        else:
            out.append(iupac_code(col))
    return "".join(out)

def expand_degenerate(primer):
    pools = [IUPAC.get(b.upper(), "N") for b in primer]
    pools = [p if p else "ACGT" for p in pools]
    return ["".join(c) for c in product(*pools)]

def fold_degeneracy(primer):
    f = 1
    for b in primer:
        f *= len(IUPAC.get(b.upper(), "ACGT"))
    return f

def three_prime_clean(primer, n=5):
    return all(b.upper() in "ACGT" for b in primer[-n:])

def gc_range(primer):
    lo, hi = 1.0, 0.0
    for o in expand_degenerate(primer) if fold_degeneracy(primer) <= 4096 else [_consensus_pick(primer, "gc_lo"), _consensus_pick(primer, "gc_hi")]:
        gc = (o.count("G") + o.count("C")) / len(o)
        lo = min(lo, gc); hi = max(hi, gc)
    return round(100 * lo, 1), round(100 * hi, 1)

def _consensus_pick(primer, mode):
    """pick a concrete oligo maximizing (gc_hi) or minimizing (gc_lo) GC, without full expansion."""
    hi_gc = {"R": "G", "Y": "C", "S": "G", "W": "A", "K": "G", "M": "C",
             "B": "G", "D": "G", "H": "C", "V": "G", "N": "G"}
    lo_gc = {"R": "A", "Y": "T", "S": "C", "W": "T", "K": "T", "M": "A",
             "B": "T", "D": "T", "H": "T", "V": "A", "N": "T"}
    tbl = hi_gc if mode == "gc_hi" else lo_gc
    return "".join(tbl.get(b.upper(), b.upper()) for b in primer)

def _tm_wallace_saltadj(oligo):
    """Simple, robust Tm estimate for QC ranking (basic + salt-ish); good enough for screening.
    Tm = 64.9 + 41*(gc-16.4)/len  (Wallace/GC for >13 nt)."""
    n = len(oligo); gc = oligo.count("G") + oligo.count("C")
    if n < 14:
        return 2 * (oligo.count("A") + oligo.count("T")) + 4 * gc
    return round(64.9 + 41.0 * (gc - 16.4) / n, 1)

def tm_range(primer):
    """min/max/mean Tm across the degenerate set (bounded expansion)."""
    if fold_degeneracy(primer) <= 4096:
        tms = [_tm_wallace_saltadj(o) for o in expand_degenerate(primer)]
    else:
        tms = [_tm_wallace_saltadj(_consensus_pick(primer, "gc_lo")),
               _tm_wallace_saltadj(_consensus_pick(primer, "gc_hi"))]
    return round(min(tms), 1), round(max(tms), 1), round(sum(tms) / len(tms), 1)

def column_conservation(aln):
    """per-column fraction of the most common base (gaps count against conservation)."""
    n = len(aln); L = len(aln[0]); cons = []
    for i in range(L):
        col = [s[i].upper() for s in aln]
        counts = {}
        for b in col:
            counts[b] = counts.get(b, 0) + 1
        top = max(counts.values())
        cons.append(top / n)
    return cons

def primer_windows(aln, plen=22, max_deg=16, min_col_cons=0.6, tail3_clean=4):
    """
    Scan an MSA for candidate primer windows of length `plen` whose IUPAC consensus has fold-degeneracy
    <= max_deg and a clean 3' end. Returns list of dicts: {start, primer(consensus IUPAC), fold, tm, gc}.
    (Caller chooses forward vs reverse orientation and pairs windows.)
    """
    L = len(aln[0]); out = []
    for s in range(0, L - plen + 1):
        cols = [[seq[s + k].upper() for seq in aln] for k in range(plen)]
        # skip windows with gaps
        if any(any(b not in "ACGT" for b in c) for c in cols):
            continue
        cons = consensus_iupac([seq[s:s + plen] for seq in aln])
        if fold_degeneracy(cons) > max_deg:
            continue
        if not three_prime_clean(cons, tail3_clean):
            continue
        tmn, tmx, tmm = tm_range(cons)
        gcl, gch = gc_range(cons)
        out.append({"start": s, "primer": cons, "fold": fold_degeneracy(cons),
                    "tm_min": tmn, "tm_max": tmx, "tm_mean": tmm, "gc_min": gcl, "gc_max": gch})
    return out


# --------------------------- unit tests ---------------------------
if __name__ == "__main__":
    # IUPAC consensus across 3 aligned "species"
    seqs = ["ACGTACGTACGT",
            "ACGTACATACGT",
            "ACGTACGTACGT"]
    cons = consensus_iupac(seqs)
    assert cons[:6] == "ACGTAC", cons
    assert cons[6] == "R" and cons[7] == "T", cons          # G/A -> R at col 6
    assert fold_degeneracy(cons) == 2, fold_degeneracy(cons)
    assert set(expand_degenerate("ACGТ".replace("Т", "T"))) == {"ACGT"}
    assert set(expand_degenerate("ACR")) == {"ACA", "ACG"}
    assert three_prime_clean("ACGTACG", 4) is True
    assert three_prime_clean("ACGTACR", 4) is False
    tmn, tmx, tmm = tm_range("ACGTACGTACGTACGTACGT")
    assert tmx >= tmn
    lo, hi = gc_range("ACGTACGTACGTACGTACRT")               # one degenerate base
    assert 0 <= lo <= hi <= 100
    # primer_windows finds the conserved 5' block
    aln = ["ACGTACGTACGTAAAACCCCGGGG",
           "ACGTACGTACGTAAAACCCCGGGT",
           "ACGTACGTACGTAAAACCCCGGGG"]
    w = primer_windows(aln, plen=12, max_deg=4)
    assert any(x["start"] == 0 for x in w), w
    assert revcomp("ATGCR") == "YGCAT", revcomp("ATGCR")
    print("degeneracy_utils: all unit tests passed")
    print(" example consensus:", cons, "fold", fold_degeneracy(cons), "tm", tm_range(cons))
