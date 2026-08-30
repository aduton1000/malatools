#!/usr/bin/env python3
"""
In-silico PCR against reference genomes — the core validation engine.

Given IUPAC primers, finds binding sites (mismatch-tolerant, with a clean 3' requirement), predicts
products across every reference sequence, and returns per-target amplicons + insert sequences. Used to
(a) validate published primers, (b) validate our own designs, (c) check specificity.
"""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from degeneracy_utils import IUPAC, revcomp

def _match_at(primer, seq, i, max_mm, clean3):
    """number of mismatches of primer vs seq starting at i, or None if > max_mm or 3'-end dirty.
    primer 3' end = right side."""
    mm = 0
    L = len(primer)
    for k in range(L):
        t = seq[i + k]
        allowed = IUPAC.get(primer[k], "")
        if t not in allowed:
            mm += 1
            if k >= L - clean3:          # mismatch within the 3' clean zone
                return None
            if mm > max_mm:
                return None
    return mm

def find_sites(primer, seq, max_mm=3, clean3=3):
    """positions where primer matches seq (+ strand) within tolerance. Returns [(pos, mm)]."""
    out = []
    L = len(primer)
    for i in range(len(seq) - L + 1):
        mm = _match_at(primer, seq, i, max_mm, clean3)
        if mm is not None:
            out.append((i, mm))
    return out

def ispcr(fwd, rev, seq, min_prod=80, max_prod=2000, max_mm=3, clean3=3):
    """Return list of products {start,end,size,fwd_mm,rev_mm,insert_start,insert_end}.
    fwd binds + strand; rev is the reverse primer as ordered (its revcomp appears on + strand)."""
    prods = []
    rev_rc = revcomp(rev)
    f_sites = find_sites(fwd, seq, max_mm, clean3)
    r_sites = find_sites(rev_rc, seq, max_mm, clean3)     # revcomp of rev on + strand
    for fp, fmm in f_sites:
        for rp, rmm in r_sites:
            end = rp + len(rev_rc)
            size = end - fp
            if min_prod <= size <= max_prod and rp >= fp + len(fwd):
                prods.append({"start": fp, "end": end, "size": size, "fwd_mm": fmm, "rev_mm": rmm,
                              "ins_start": fp + len(fwd), "ins_end": rp})
    return prods

def ispcr_bothstrands(fwd, rev, seq, **kw):
    """try the amplicon on + strand and on the reverse complement of the template."""
    res = [("+", p) for p in ispcr(fwd, rev, seq, **kw)]
    rc = revcomp(seq)
    res += [("-", p) for p in ispcr(fwd, rev, rc, **kw)]
    return res

def load_fasta(path):
    d = {}
    for r in open(path).read().split(">"):
        if r.strip():
            h, s = r.split("\n", 1); d[h.split()[0]] = s.replace("\n", "").upper()
    return d

def validate_pair(name, fwd, rev, refs, order, max_mm=4, clean3=3, min_prod=80, max_prod=2000):
    """Run isPCR of (fwd,rev) against each reference; print a report + return per-target inserts."""
    print(f"\n=== {name} ===\n  FWD {fwd}\n  REV {rev}")
    inserts = {}
    for sp in order:
        if sp not in refs:
            print(f"  {sp:>4}: (no ref)"); continue
        hits = ispcr_bothstrands(fwd, rev, refs[sp], max_mm=max_mm, clean3=clean3,
                                 min_prod=min_prod, max_prod=max_prod)
        if not hits:
            print(f"  {sp:>4}: NO PRODUCT"); continue
        strand, p = min(hits, key=lambda x: x[1]["size"])         # smallest product = the intended one
        s = refs[sp] if strand == "+" else revcomp(refs[sp])
        insert = s[p["ins_start"]:p["ins_end"]]
        inserts[sp] = insert
        print(f"  {sp:>4}: {p['size']:>4} bp  (fwd_mm {p['fwd_mm']}, rev_mm {p['rev_mm']}, strand {strand}, insert {len(insert)} bp)")
    return inserts

def discrimination(inserts, order):
    """align-free pairwise SNP counts on equal-length inserts (fallback: length diff)."""
    keys = [k for k in order if k in inserts]
    print("  pairwise insert distance (SNP if equal length, else 'L±'):")
    print("       " + "  ".join(f"{k:>4}" for k in keys))
    for i, a in enumerate(keys):
        row = []
        for j, b in enumerate(keys):
            ia, ib = inserts[a], inserts[b]
            if len(ia) == len(ib):
                row.append(f"{sum(1 for x,y in zip(ia,ib) if x!=y):>4}")
            else:
                row.append(f"{'L'+('+' if len(ia)>len(ib) else '-'):>4}")
        print(f"  {a:>4} " + "  ".join(row))
    mn = min((sum(1 for x,y in zip(inserts[a],inserts[b]) if x!=y) if len(inserts[a])==len(inserts[b]) else 999)
             for i,a in enumerate(keys) for b in keys[i+1:]) if len(keys)>1 else 0
    return mn

if __name__ == "__main__":
    refs = load_fasta("refs/mito/pan_mito.labeled.fasta")
    order = ["Pf", "Pv", "Poc", "Pow", "Pm", "Pk"]
    # published primers to validate
    tests = [
        ("cox3 Echeverry short",  "AGCGGTTAACCTTTCYTTTTCCTTACG", "AGTGCATCATGTATGACAGCATGTTTACA"),
        ("cox3 Isozumi MtUF/MtUR","CTCGCCATTTGATAGCGGTTAACC",    "CCTGTTATCCCCGGCGAACCTTC"),
        ("cytb Xu 430",           "TGGTAGCACAAATCCTTTAGGG",      "TGGTAATTGACATCCAATCC"),
        ("cox1 Echeverry 2017",   "AGAACGAACGCTTTTAACGCCTG",     "ACTTAATGGTGGATATAAAGTCCATCCWGT"),
    ]
    for name, f, r in tests:
        ins = validate_pair(name, f, r, refs, order, max_mm=5, clean3=3, min_prod=120, max_prod=1200)
        if len(ins) >= 2:
            discrimination(ins, order)
