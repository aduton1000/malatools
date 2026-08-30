#!/usr/bin/env python3
"""
Locate resistance-gene orthologs in each species genome by tblastn of the Pf reference protein.
Robust to product-name variation and to intron structure (HSPs are stitched into a gene span).
Outputs targets/orthologs.tsv: gene, species, contig, start, end, strand, pct_id, hsp_cov.
"""
import os
import gzip, subprocess, re, os
import shutil
# Bin dir for BLAST+/mafft — resolved from $MALAFORGE_BIN or PATH (no hardcoded machine path).
BB = (os.environ.get("MALAFORGE_BIN")
      or (os.path.dirname(shutil.which("blastn") or shutil.which("makeblastdb")
                          or shutil.which("mafft") or shutil.which("tblastn") or "") or ""))
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (design/ -> ..)

# Pf reference resistance genes (locus tags in Pf3D7.gff)
PF_GENES = {
    "dhfr": "PF3D7_0417200", "dhps": "PF3D7_0810800", "crt": "PF3D7_0709000",
    "mdr1": "PF3D7_0523000", "k13": "PF3D7_1343700",
}
SPECIES_GENOME = {
    "Pm": "refs/nuclear/PmUG01.fna", "Pk": "refs/nuclear/PkHl.fna",
    "Poc": "refs/nuclear/PocGH01.fna", "Pow": "refs/nuclear/PowCR01.fna",
    "Pv": "refs/nuclear/PvP01.fna",
}
bases = "TCAG"; aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
C2A = {c1 + c2 + c3: aas[i * 16 + j * 4 + k] for i, c1 in enumerate(bases)
       for j, c2 in enumerate(bases) for k, c3 in enumerate(bases)}
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
def translate(nt):
    return "".join(C2A.get(nt[i:i+3], "X") for i in range(0, len(nt) - 2, 3))

def load_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    d = {}; name = None; buf = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name: d[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if name: d[name] = "".join(buf)
    return d

def attr(s, k):
    m = re.search(rf'{k}=([^;]+)', s); return m.group(1) if m else None

def extract_pf_proteins():
    """parse Pf3D7.gff (mRNA->CDS) and translate the 5 genes."""
    fa = load_fasta(f"{BASE}/refs/nuclear/Pf3D7.fna")
    want = set(PF_GENES.values()); rid2lt = {}; cds = {}
    with gzip.open(f"{BASE}/refs/anno/Pf3D7.gff.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            if f[2] == "mRNA":
                par = attr(f[8], "Parent"); rid = attr(f[8], "ID")
                if par and par.startswith("gene-") and par[5:] in want and par[5:] not in rid2lt.values():
                    rid2lt[rid] = par[5:]
            elif f[2] == "CDS":
                par = attr(f[8], "Parent")
                if par in rid2lt:
                    cds.setdefault(rid2lt[par], {"contig": f[0], "strand": f[6], "segs": []})
                    cds[rid2lt[par]]["segs"].append((int(f[3]), int(f[4])))
    prots = {}
    lt2gene = {v: k for k, v in PF_GENES.items()}
    for lt, g in cds.items():
        segs = sorted(g["segs"], reverse=(g["strand"] == "-"))
        nt = ""
        for s, e in segs:
            piece = fa[g["contig"]][s - 1:e]
            nt += rc(piece) if g["strand"] == "-" else piece
        prots[lt2gene[lt]] = translate(nt).rstrip("*")
    return prots

def main():
    os.chdir(BASE)
    prots = extract_pf_proteins()
    print("Pf reference proteins:", {k: len(v) for k, v in prots.items()})
    # build blastdbs
    for sp, gpath in SPECIES_GENOME.items():
        subprocess.run([f"{BB}/makeblastdb", "-in", gpath, "-dbtype", "nucl",
                        "-out", f"scratch/db_{sp}"], capture_output=True)
    rows = [("gene", "species", "contig", "start", "end", "strand", "pct_id", "aa_cov")]
    for gene, prot in prots.items():
        open("scratch/q.faa", "w").write(f">{gene}\n{prot}\n")
        for sp in SPECIES_GENOME:
            r = subprocess.run([f"{BB}/tblastn", "-query", "scratch/q.faa", "-db", f"scratch/db_{sp}",
                                "-outfmt", "6 sseqid sstart send sframe pident length", "-max_target_seqs", "5"],
                               capture_output=True, text=True)
            hsps = [l.split("\t") for l in r.stdout.strip().splitlines()]
            if not hsps: rows.append((gene, sp, "NONE", "", "", "", "", "")); continue
            # anchor on the best HSP (first row) and keep only HSPs of the SAME locus
            # (same contig AND within ~10 kb of the anchor midpoint) -> excludes distant paralogs
            top = hsps[0]; tc = top[0]; tmid = (int(top[1]) + int(top[2])) / 2
            same = [h for h in hsps if h[0] == tc and abs((int(h[1]) + int(h[2])) / 2 - tmid) < 10000]
            starts = [int(h[1]) for h in same]; ends = [int(h[2]) for h in same]
            allpos = starts + ends
            frame = int(top[3]); strand = "+" if frame > 0 else "-"
            pid = max(float(h[4]) for h in same); cov = sum(int(h[5]) for h in same)
            rows.append((gene, sp, tc, min(allpos), max(allpos), strand, f"{pid:.1f}", f"{cov}/{len(prot)}"))
    with open("targets/orthologs.tsv", "w") as fh:
        for r in rows: fh.write("\t".join(map(str, r)) + "\n")
    for r in rows: print("  " + "\t".join(map(str, r)))

if __name__ == "__main__":
    main()
