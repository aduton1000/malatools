#!/usr/bin/env python3
"""
panplasmo.py — end-to-end analysis driver for one sample (map -> all analysis stages).
Mirrors the Nextflow workflow; used for local runs and CI validation.

Usage: panplasmo.py --fastq R.fastq.gz --sample NAME --outdir OUT \
         [--ref pan_reference.fa.gz] [--manifest amplicons.tsv] [--codon-map ...] [--pk-panel ...]
         [--platform ont|illumina] [--r2 R2.fastq.gz]
"""
import sys, os, argparse, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, "../.."))
PY = sys.executable


def sh(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def run_capture_species(callfile):
    sp = open(callfile).read().strip().split("\t") if os.path.exists(callfile) else []
    return [s for s in sp if s]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True); ap.add_argument("--r2")
    ap.add_argument("--sample", required=True); ap.add_argument("--outdir", required=True)
    ap.add_argument("--ref", default=f"{BASE}/refs/pan_reference.fa.gz")
    ap.add_argument("--manifest", default=f"{BASE}/analysis/manifest/amplicons.tsv")
    ap.add_argument("--codon-map", default=f"{BASE}/targets/codon_to_genomic_pan.tsv")
    ap.add_argument("--pk-panel", default=f"{BASE}/analysis/refs/pk_diag_panel.tsv")
    ap.add_argument("--platform", default="ont", choices=["ont", "illumina"])
    ap.add_argument("--threads", default="4")
    a = ap.parse_args()
    od = os.path.join(a.outdir, a.sample); os.makedirs(od, exist_ok=True)
    bam = f"{od}/{a.sample}.bam"
    env = dict(os.environ, PYTHONPATH=HERE)

    # 1. MAP (minimap2 -> sorted, indexed BAM)
    preset = "map-ont" if a.platform == "ont" else "sr"
    mm = ["minimap2", "-ax", preset, "-t", a.threads, "--secondary=no", a.ref, a.fastq]
    if a.r2:
        mm.append(a.r2)
    p1 = subprocess.Popen(mm, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["samtools", "sort", "-@", a.threads, "-o", bam], stdin=p1.stdout, stderr=subprocess.DEVNULL)
    p1.stdout.close(); p2.communicate()
    sh(["samtools", "index", bam])

    # 2. SPECIES ID
    sp_tsv = f"{od}/{a.sample}.species.tsv"
    sh([PY, f"{HERE}/species_id.py", bam, a.manifest, a.sample, sp_tsv], env=env)
    present = run_capture_species(sp_tsv.replace(".tsv", ".call"))
    sp_arg = ",".join(present) if present else "Pf,Pv,Poc,Pow,Pm,Pk"

    # 3. RESISTANCE
    sh([PY, f"{HERE}/resistance_genotype.py", bam, a.manifest, a.codon_map, a.sample, od, "--species", sp_arg], env=env)
    # 4. CNV / DELETIONS
    sh([PY, f"{HERE}/cnv_deletion.py", bam, a.manifest, a.sample, f"{od}/{a.sample}.cnv.tsv", "--species", sp_arg], env=env)
    # 5. DIVERSITY
    sh([PY, f"{HERE}/diversity.py", bam, a.manifest, a.sample, od, "--species", sp_arg], env=env)
    # 6. PK CLUSTER (only if Pk present)
    if "Pk" in present and os.path.exists(a.pk_panel):
        sh([PY, f"{HERE}/pk_cluster.py", bam, a.manifest, a.pk_panel, a.sample, f"{od}/{a.sample}.cluster.tsv"], env=env)
    # 7. REPORT
    if os.path.exists(f"{HERE}/report.py"):
        sh([PY, f"{HERE}/report.py", od, a.sample], env=env)
    print(f"[panplasmo] {a.sample}: done -> {od}  (species: {sp_arg})")


if __name__ == "__main__":
    main()
