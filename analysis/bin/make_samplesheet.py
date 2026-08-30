#!/usr/bin/env python3
"""
make_samplesheet.py — turn a platform-native sequencing run folder into a ready-to-use samplesheet
(`sample,fastq[,fastq2]`) for the pan-plasmodium pipeline.

Understands the three run-folder layouts:
  ONT (MinKNOW)   run/fastq_pass/barcodeNN/*.fastq.gz   (one barcode = one sample; also 'unclassified')
                  or a single-sample run: run/fastq_pass/*.fastq.gz
  Illumina        run/**/SampleName_S1_L001_R1_001.fastq.gz  (+ _R2_)   [bcl2fastq / BCL Convert]
  MGI (DNBSEQ)    run/**/SampleName_1.fq.gz  (+ _2.fq.gz)   or  *_R1.fastq.gz / *_R2.fastq.gz

Multiple files per sample (ONT chunks, Illumina/MGI lanes) are concatenated into one file per sample
(gzip concatenation is stream-safe — no re-compression) in --merge-dir, so the emitted paths are single,
ready-to-consume files. Use --barcode-map to rename ONT barcodes to real sample IDs.

Examples:
  make_samplesheet.py --platform ont      --run-dir /data/minknow/RUN1        --out samples.csv
  make_samplesheet.py --platform illumina --run-dir /data/miseq/240830_M0/Fastq --out samples.csv
  make_samplesheet.py --platform mgi      --run-dir /data/mgi/E100xxx          --out samples.csv
  make_samplesheet.py --platform ont --run-dir RUN1 --barcode-map barcodes.csv --out samples.csv
"""
import os, sys, re, csv, glob, argparse, subprocess

FASTQ_EXT = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def find_fastqs(root):
    out = []
    for dp, _, fns in os.walk(root):
        for f in fns:
            if f.endswith(FASTQ_EXT) and not f.startswith("."):
                out.append(os.path.join(dp, f))
    return sorted(out)


def merge(files, dest):
    """concatenate 1+ fastq(.gz) into dest (gzip-concat safe); return dest. Single file -> just abspath."""
    files = sorted(files)
    if len(files) == 1:
        return os.path.abspath(files[0])
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "wb") as out:
        for f in files:
            with open(f, "rb") as fh:
                while True:
                    chunk = fh.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
    return os.path.abspath(dest)


def load_barcode_map(path):
    m = {}
    if path and os.path.exists(path):
        for r in csv.reader(open(path)):
            if len(r) >= 2 and r[0].strip().lower() not in ("barcode", ""):
                m[r[0].strip()] = r[1].strip()
    return m


# ---------- ONT ----------
def parse_ont(run, mdir, bmap, min_bytes, keep_unclassified):
    # locate the reads root (fastq_pass if present, else the run dir)
    roots = glob.glob(os.path.join(run, "**", "fastq_pass"), recursive=True) or [run]
    root = roots[0]
    samples = {}
    bc_dirs = sorted(d for d in glob.glob(os.path.join(root, "*")) if os.path.isdir(d)
                     and re.match(r"(barcode\d+|unclassified|BC\d+)$", os.path.basename(d), re.I))
    if bc_dirs:
        for d in bc_dirs:
            bc = os.path.basename(d)
            if bc.lower() == "unclassified" and not keep_unclassified:
                continue
            fq = [f for f in find_fastqs(d)]
            if sum(os.path.getsize(f) for f in fq) < min_bytes:
                sys.stderr.write(f"  skip {bc}: empty/below --min-bytes\n"); continue
            sample = bmap.get(bc, bc)
            samples[sample] = merge(fq, os.path.join(mdir, f"{sample}.fastq.gz"))
    else:  # single-sample run
        fq = find_fastqs(root)
        if fq:
            sample = bmap.get("sample", os.path.basename(os.path.normpath(run)))
            samples[sample] = merge(fq, os.path.join(mdir, f"{sample}.fastq.gz"))
    return samples, None


# ---------- Illumina / MGI (paired) ----------
R1_TOKENS = [r"_R1_001", r"_R1(?=\.)", r"_1(?=\.f(ast)?q)", r"_R1(?=_)"]
def r1_to_r2(name):
    for a, b in [("_R1_001", "_R2_001"), ("_R1.", "_R2."), ("_R1_", "_R2_"), ("_1.fq", "_2.fq"),
                 ("_1.fastq", "_2.fastq")]:
        if a in name:
            return name.replace(a, b, 1)
    return None


def sample_name(fname):
    b = re.sub(r"\.(fastq|fq)(\.gz)?$", "", os.path.basename(fname))
    b = re.sub(r"_R[12]_001$", "", b)
    b = re.sub(r"_R[12]$", "", b)
    b = re.sub(r"_[12]$", "", b)
    b = re.sub(r"_S\d+(_L\d{3})?$", "", b)         # Illumina _S1 / _S1_L001
    b = re.sub(r"_L\d{2,3}$", "", b)               # lane
    b = re.sub(r"_S\d+_L\d{3}", "", b)
    return b or os.path.basename(fname)


def parse_paired(run, mdir, min_bytes):
    all_fq = find_fastqs(run)
    is_r1 = lambda f: bool(re.search(r"(_R1_001|_R1\.|_R1_|_1\.f(ast)?q)", os.path.basename(f)))
    r1s = [f for f in all_fq if is_r1(f)]
    by_sample = {}
    unpaired = []
    for r1 in r1s:
        r2 = r1_to_r2(r1)
        r2 = r2 if (r2 and os.path.exists(r2)) else None
        s = sample_name(r1)
        by_sample.setdefault(s, {"r1": [], "r2": []})
        by_sample[s]["r1"].append(r1)
        if r2:
            by_sample[s]["r2"].append(r2)
        else:
            unpaired.append(os.path.basename(r1))
    samples = {}
    for s, d in by_sample.items():
        if sum(os.path.getsize(f) for f in d["r1"]) < min_bytes:
            sys.stderr.write(f"  skip {s}: empty/below --min-bytes\n"); continue
        r1 = merge(d["r1"], os.path.join(mdir, f"{s}.R1.fastq.gz"))
        r2 = merge(d["r2"], os.path.join(mdir, f"{s}.R2.fastq.gz")) if d["r2"] else ""
        samples[s] = (r1, r2)
    return samples, unpaired


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--platform", required=True, choices=["ont", "illumina", "mgi"])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--merge-dir", default=None, help="where to write merged fastqs (default: <out_dir>/merged)")
    ap.add_argument("--barcode-map", default=None, help="ONT: CSV barcode,sample to rename barcodes")
    ap.add_argument("--min-bytes", type=int, default=1000, help="skip samples smaller than this (empty barcodes)")
    ap.add_argument("--keep-unclassified", action="store_true")
    a = ap.parse_args()
    if not os.path.isdir(a.run_dir):
        sys.exit(f"run-dir not found: {a.run_dir}")
    mdir = a.merge_dir or os.path.join(os.path.dirname(os.path.abspath(a.out)) or ".", "merged")

    if a.platform == "ont":
        samples, unpaired = parse_ont(a.run_dir, mdir, load_barcode_map(a.barcode_map), a.min_bytes, a.keep_unclassified)
        header = ["sample", "fastq"]
        rows = [[s, samples[s]] for s in sorted(samples)]
    else:
        samples, unpaired = parse_paired(a.run_dir, mdir, a.min_bytes)
        header = ["sample", "fastq", "fastq2"]
        rows = [[s, samples[s][0], samples[s][1]] for s in sorted(samples)]
        if unpaired:
            sys.stderr.write(f"  WARNING: {len(unpaired)} R1 file(s) had no R2 mate (treated single-end): "
                             + ", ".join(unpaired[:4]) + ("..." if len(unpaired) > 4 else "") + "\n")

    if not rows:
        sys.exit(f"No FASTQ samples found under {a.run_dir} for platform {a.platform}. "
                 "Check the path / layout, or --keep-unclassified for ONT.")
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    sys.stderr.write(f"[make_samplesheet] {a.platform}: {len(rows)} samples -> {a.out}"
                     + (f"  (merged fastqs in {mdir})" if os.path.isdir(mdir) else "") + "\n")
    for r in rows[:8]:
        sys.stderr.write("   " + " ".join(os.path.basename(x) for x in r) + "\n")
    print(a.out)


if __name__ == "__main__":
    main()
