#!/usr/bin/env python3
"""malacall — analysis & reporting pipeline for the pan-Plasmodium amplicon panel.

One blood-spot in → species ID (6 taxa, mixed infections), drug-resistance genotypes, RDT-escape
deletions & gene amplification, diversity/COI, P. knowlesi subpopulation, per-sample clinical reports,
and a cohort surveillance roll-up with native COI-MLE + relatedness.
"""
import os, sys, argparse, subprocess
from ._common import (REPO, run_script, find_nextflow, java_home, dispatch, _c, C_TEAL, C_DIM)

TAGLINE = "pan-Plasmodium analysis & reporting pipeline"
EPILOG = ("MalaCall wraps the Nextflow workflow and the per-stage Python tools. Pair it with MalaForge, "
          "the panel-design toolkit. Docs: analysis/MANUAL.md")

COMMANDS = [
    ("run",         "Run the full pipeline on a samplesheet (Nextflow; the usual entry point)"),
    ("sample",      "Run ONE sample end-to-end without Nextflow (pure Python driver)"),
    ("samplesheet", "Build a samplesheet from a platform run folder (ONT/Illumina/MGI)"),
    ("demux",       "Demultiplex a pooled FASTQ by the panel's custom dual barcodes"),
    ("manifest",    "Build the amplicon manifest on the pan-reference (one-time setup)"),
    ("cohort",      "Cohort surveillance roll-up across a results directory"),
    ("relatedness", "Cohort COI (MLE) + pairwise relatedness (+ McCOIL/Dcifer exports)"),
    ("report",      "(Re)render one sample's clinical HTML report"),
    ("validate",    "Run the simulated-truth end-to-end validation (62 assertions)"),
]


def cmd_run(argv):
    """Smart wrapper over `nextflow run analysis/main.nf …` — sets JAVA_HOME, defaults the profile,
    and translates friendly --dashed flags to Nextflow's underscore params."""
    p = argparse.ArgumentParser(
        prog="malacall run", description="Run the full pipeline on a samplesheet (via Nextflow).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--samplesheet", required=True, help="samplesheet CSV (see `malacall samplesheet`)")
    p.add_argument("--outdir", default="results", help="output directory")
    p.add_argument("--platform", default="ont", choices=["ont", "illumina", "mgi"], help="sequencing platform")
    p.add_argument("--profile", default="localpy", help="Nextflow profile (env with pysam/numpy/pandas)")
    p.add_argument("--cohort-name", help="name for the cohort roll-up outputs")
    p.add_argument("--min-frac", type=float, help="species-call read fraction")
    p.add_argument("--min-vaf", type=float, help="minor resistant-clone VAF")
    p.add_argument("--resume", action="store_true", help="resume a previous run")
    # demultiplexing
    g = p.add_argument_group("demultiplexing (optional)")
    g.add_argument("--demux", action="store_true", help="demultiplex a pooled FASTQ first")
    g.add_argument("--pooled", help="glob of pooled FASTQs, e.g. 'run/fastq_pass/*.fastq.gz'")
    g.add_argument("--barcode-map", help="sample,fwd_bc,rev_bc CSV")
    g.add_argument("--barcodes", help="barcode sequences TSV (default: bundled 96-set)")
    g.add_argument("--demux-mode", choices=["dual", "single"], help="barcode scheme")
    g.add_argument("--bc-mismatch", type=int, help="barcode mismatch tolerance")
    # advanced / plumbing
    a = p.add_argument_group("advanced")
    a.add_argument("--manifest", help="amplicon manifest TSV")
    a.add_argument("--ref", help="pan-reference FASTA(.gz)")
    a.add_argument("--codon-map", help="codon→genomic map TSV")
    a.add_argument("--pk-panel", help="Pk cluster diagnostic panel TSV")
    a.add_argument("--python", help="interpreter Nextflow runs the stages with (default: this CLI's own)")
    a.add_argument("--work-dir", help="Nextflow work directory")
    a.add_argument("--nextflow", help="path to the nextflow launcher")
    a.add_argument("--java-home", help="JAVA_HOME to use for Nextflow")
    a.add_argument("--dry-run", action="store_true", help="print the nextflow command and exit")
    a.add_argument("extra", nargs=argparse.REMAINDER,
                   help="everything after `--` is passed straight to nextflow")
    args = p.parse_args(argv)

    nf = find_nextflow(args.nextflow)
    if not nf:
        sys.exit("malacall run: nextflow not found. Install it, put it on PATH, pass --nextflow PATH, "
                 "or set MALACALL_NEXTFLOW. (Or use `malacall sample` for a Nextflow-free single-sample run.)")
    main_nf = os.path.join(REPO, "analysis", "main.nf")

    # friendly flag -> nextflow --underscore_param
    params = {
        "samplesheet": args.samplesheet, "outdir": args.outdir, "platform": args.platform,
        "cohort_name": args.cohort_name, "min_frac": args.min_frac, "min_vaf": args.min_vaf,
        "pooled": args.pooled, "barcode_map": args.barcode_map, "barcodes": args.barcodes,
        "demux_mode": args.demux_mode, "bc_mismatch": args.bc_mismatch,
        "manifest": args.manifest, "ref": args.ref, "codon_map": args.codon_map, "pk_panel": args.pk_panel,
        # inject the exact interpreter so nextflow.config needs no hardcoded machine path
        "python": args.python or sys.executable,
    }
    cmd = [nf, "run", main_nf, "-profile", args.profile]
    if args.resume:
        cmd.append("-resume")
    if args.work_dir:
        cmd += ["-work-dir", args.work_dir]
    if args.demux:
        cmd.append("--demux")
    for k, v in params.items():
        if v is not None:
            cmd += [f"--{k}", str(v)]
    extra = args.extra[1:] if args.extra and args.extra[0] == "--" else args.extra
    if extra:
        cmd += extra

    env = dict(os.environ)
    jh = java_home(args.java_home)
    if jh:
        env["JAVA_HOME"] = jh
        env["JAVA_CMD"] = os.path.join(jh, "bin", "java")

    printable = " ".join(cmd)
    print(_c("malacall run →", C_TEAL), printable)
    if jh:
        print(_c(f"  (JAVA_HOME={jh})", C_DIM))
    if args.dry_run:
        return 0
    return subprocess.call(cmd, cwd=os.path.join(REPO, "analysis"), env=env)


def cmd_validate(argv):
    """Regenerate simulated truth, run the pipeline on every simulated sample, and score it."""
    p = argparse.ArgumentParser(prog="malacall validate",
                                description="End-to-end validation on bundled simulated data.")
    p.add_argument("--sim", default=os.path.join(REPO, "analysis", "test", "sim"), help="sim working dir")
    p.add_argument("--out", default=os.path.join(REPO, "analysis", "test", "out"), help="results dir")
    p.add_argument("--skip-regen", action="store_true", help="score existing --out without re-running")
    a = p.parse_args(argv)
    py = sys.executable
    A = os.path.join(REPO, "analysis")
    ref = os.path.join(REPO, "refs", "pan_reference.fa.gz")
    man = os.path.join(A, "manifest", "amplicons.tsv")
    codon = os.path.join(REPO, "targets", "codon_to_genomic_pan.tsv")
    pk = os.path.join(A, "refs", "pk_diag_panel.tsv")
    if not a.skip_regen:
        os.makedirs(a.sim, exist_ok=True)
        print(_c("[validate] simulating known-truth reads…", C_DIM))
        if subprocess.call([py, os.path.join(A, "test", "simulate.py"), man, ref, codon, pk, a.sim]):
            return 1
        reads = os.path.join(a.sim, "reads")
        for fq in sorted(os.listdir(reads)):
            if not fq.endswith(".fastq.gz"):
                continue
            s = fq[:-len(".fastq.gz")]
            print(_c(f"[validate] pipeline → {s}", C_DIM))
            subprocess.call([py, os.path.join(A, "bin", "panplasmo.py"),
                             "--fastq", os.path.join(reads, fq), "--sample", s,
                             "--outdir", a.out, "--platform", "ont"])
        subprocess.call([py, os.path.join(A, "bin", "cohort_report.py"), a.out, "validation"])
        subprocess.call([py, os.path.join(A, "bin", "relatedness.py"), a.out, "validation"])
    return subprocess.call([py, os.path.join(A, "test", "validate.py"), a.sim, a.out])


HANDLERS = {
    "run": cmd_run,
    "sample": "analysis/bin/panplasmo.py",
    "samplesheet": "analysis/bin/make_samplesheet.py",
    "demux": "analysis/bin/demux.py",
    "manifest": "analysis/bin/build_manifest.py",
    "cohort": "analysis/bin/cohort_report.py",
    "relatedness": "analysis/bin/relatedness.py",
    "report": "analysis/bin/report.py",
    "validate": cmd_validate,
}


def main():
    sys.exit(dispatch("malacall", TAGLINE, COMMANDS, HANDLERS, EPILOG))


if __name__ == "__main__":
    main()
