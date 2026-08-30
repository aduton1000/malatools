#!/usr/bin/env python3
"""
species_id.py — detect & speciate all 6 taxa from mito cox3 read proportions; flag mixed infections.
Usage: species_id.py <bam> <manifest.tsv> <sample> <out.tsv> [--min-reads N] [--min-frac F]
Output: per-species read count, proportion, and call (present/absent), + a species summary line.
"""
import sys, csv, argparse
import pysam
from common import load_manifest, region_read_count

ALL = ["Pf", "Pv", "Poc", "Pow", "Pm", "Pk"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam"); ap.add_argument("manifest"); ap.add_argument("sample"); ap.add_argument("out")
    ap.add_argument("--min-reads", type=int, default=20)
    ap.add_argument("--min-frac", type=float, default=0.01)
    a = ap.parse_args()
    man = load_manifest(a.manifest)
    cox3 = [m for m in man if m["role"] == "species_id"]
    bam = pysam.AlignmentFile(a.bam, "rb")
    reads = {sp: 0 for sp in ALL}
    # cox3 reads are only ~8 SNP apart across taxa -> minimap2 gives MAPQ 0 (near-ties); they are
    # species-diagnostic by WHICH mito contig they map to, so count primary alignments at min_mapq 0.
    for m in cox3:
        sp = m["species"]
        reads[sp] += region_read_count(bam, m["pan_contig"], m["start"], m["end"], min_mapq=0)
    total = sum(reads.values())
    rows = []
    present = []
    for sp in ALL:
        frac = reads[sp] / total if total else 0.0
        called = reads[sp] >= a.min_reads and frac >= a.min_frac
        if called:
            present.append((sp, frac))
        rows.append((sp, reads[sp], round(frac, 4), "present" if called else "absent"))
    present.sort(key=lambda x: -x[1])
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "cox3_reads", "proportion", "call"])
        for sp, n, fr, c in rows:
            w.writerow([a.sample, sp, n, fr, c])
    summary = ("MIXED: " if len(present) > 1 else "") + ", ".join(
        f"{sp} ({fr*100:.0f}%)" for sp, fr in present) if present else "NO PARASITE DETECTED"
    sys.stderr.write(f"[species_id] {a.sample}: {summary}  (total cox3 reads={total})\n")
    # machine-readable call file
    with open(a.out.replace(".tsv", ".call"), "w") as fh:
        fh.write("\t".join(sp for sp, _ in present) + "\n")


if __name__ == "__main__":
    main()
