#!/usr/bin/env python3
"""
fetch_genomes.py — download the raw per-species reference genomes (+ annotations) that MalaForge's
*design* step needs. These are git-ignored because they are large (~250 MB) and public; the analysis
pipeline (MalaCall) only needs the small committed `refs/pan_reference.fa.gz`, so you only need these
if you want to re-run or extend the panel design.

Source: NCBI Datasets API (no `datasets` CLI required — plain HTTPS + zip). Files land as
`refs/nuclear/<name>.fna` and `refs/anno/<name>.gff`, matching the names the design scripts expect.

Usage:
  python3 refs/fetch_genomes.py            # download any that are missing
  python3 refs/fetch_genomes.py --force    # re-download even if present
  python3 refs/fetch_genomes.py --list     # show the plan and exit (no network)
  python3 refs/fetch_genomes.py --only Pf3D7,PkHl

VERIFY: the accessions below reproduce the assemblies used to build this panel. Confirm them against
your refs provenance / NCBI before relying on exact coordinates for a fresh design run.
"""
import os, sys, io, zipfile, argparse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))       # .../malaforge/refs
NUC  = os.path.join(HERE, "nuclear")
ANNO = os.path.join(HERE, "anno")
API  = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"

# name, NCBI assembly accession, species tag, fetch GFF too?
GENOMES = [
    ("Pf3D7",   "GCF_000002765.6", "P. falciparum 3D7",          True),
    ("PvP01",   "GCA_900093555.1", "P. vivax P01",               True),
    ("PmUG01",  "GCF_900090045.1", "P. malariae UG01",           True),
    ("PkHl",    "GCF_000006355.2", "P. knowlesi strain H",       True),
    ("PocGH01", "GCA_900090035.2", "P. ovale curtisi GH01",      True),
    ("PowCR01", "GCA_900090025.2", "P. ovale wallikeri CR01",    True),
]


def plan(rows):
    print(f"{'name':9} {'accession':18} {'species':28} targets")
    for name, acc, sp, gff in rows:
        tgt = f"nuclear/{name}.fna" + (f" + anno/{name}.gff" if gff else "")
        print(f"{name:9} {acc:18} {sp:28} {tgt}")


def fetch_one(name, acc, want_gff, force):
    fna = os.path.join(NUC, f"{name}.fna")
    gff = os.path.join(ANNO, f"{name}.gff")
    if not force and os.path.exists(fna) and (not want_gff or os.path.exists(gff)):
        print(f"[skip] {name}: already present"); return True
    types = "GENOME_FASTA" + (",GENOME_GFF" if want_gff else "")
    url = f"{API.format(acc=acc)}?include_annotation_type={types}"
    print(f"[get ] {name} <- {acc} …", flush=True)
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/zip"})
        with urllib.request.urlopen(req, timeout=120) as r:
            blob = r.read()
    except Exception as e:
        print(f"[FAIL] {name}: download error: {e}"); return False
    try:
        z = zipfile.ZipFile(io.BytesIO(blob))
    except zipfile.BadZipFile:
        print(f"[FAIL] {name}: response was not a zip (check the accession / API)"); return False
    got_fna = got_gff = False
    for member in z.namelist():
        low = member.lower()
        if low.endswith((".fna", ".fasta")) and f"/{acc}/" in member:
            os.makedirs(NUC, exist_ok=True)
            with z.open(member) as src, open(fna, "wb") as dst:
                dst.write(src.read())
            got_fna = True
        elif low.endswith(".gff") and want_gff:
            os.makedirs(ANNO, exist_ok=True)
            with z.open(member) as src, open(gff, "wb") as dst:
                dst.write(src.read())
            got_gff = True
    if not got_fna:
        print(f"[FAIL] {name}: no FASTA found in archive"); return False
    msg = f"[ok  ] {name}: nuclear/{name}.fna ({os.path.getsize(fna)//1_000_000} MB)"
    if want_gff:
        msg += f" + anno/{name}.gff" if got_gff else " (WARNING: no GFF in archive)"
    print(msg)
    return True


def main():
    ap = argparse.ArgumentParser(description="Fetch raw reference genomes for MalaForge design.")
    ap.add_argument("--force", action="store_true", help="re-download even if the file exists")
    ap.add_argument("--list", action="store_true", help="print the plan and exit (no network)")
    ap.add_argument("--only", help="comma-separated subset of names (e.g. Pf3D7,PkHl)")
    a = ap.parse_args()
    rows = GENOMES
    if a.only:
        want = {x.strip() for x in a.only.split(",")}
        rows = [g for g in GENOMES if g[0] in want]
        if not rows:
            sys.exit(f"--only matched nothing; valid names: {', '.join(g[0] for g in GENOMES)}")
    if a.list:
        plan(rows); return
    os.makedirs(NUC, exist_ok=True); os.makedirs(ANNO, exist_ok=True)
    ok = sum(fetch_one(n, acc, gff, a.force) for n, acc, _, gff in rows)
    print(f"\n{ok}/{len(rows)} genomes ready in refs/nuclear (+ refs/anno).")
    if ok:
        open(os.path.join(NUC, ".download.log"), "w").write("NUCLEAR_DOWNLOADS_DONE\n")
    sys.exit(0 if ok == len(rows) else 1)


if __name__ == "__main__":
    main()
