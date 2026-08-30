#!/usr/bin/env python3
"""
cnv_deletion.py — copy-number (gene amplification) and gene-deletion (pfhrp2/pfhrp3 RDT escape) calls by
read-depth relative to same-species single-copy control amplicons.

Usage: cnv_deletion.py <bam> <manifest.tsv> <sample> <out.tsv> --species Pf,Pv \
         [--amp-thresh 1.6] [--del-thresh 0.25] [--min-control-depth 20]
"""
import sys, csv, argparse, statistics
import pysam
from common import load_manifest, region_depth

CNV_GENES = {"pfmdr1": "chloroquine/lumefantrine (amplification)", "pfpm2": "piperaquine (plasmepsin-2 amp.)",
             "pfpm3": "piperaquine (plasmepsin-3 amp.)", "pfgch1": "antifolate (GCH1 amp.)",
             "pvmdr1": "chloroquine (Pv mdr1 amp.)"}
DEL_GENES = {"pfhrp2": "HRP2-based RDT escape", "pfhrp3": "HRP3-based RDT escape"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam"); ap.add_argument("manifest"); ap.add_argument("sample"); ap.add_argument("out")
    ap.add_argument("--species", default="Pf,Pv,Poc,Pow,Pm,Pk")
    ap.add_argument("--amp-thresh", type=float, default=1.6)
    ap.add_argument("--del-thresh", type=float, default=0.25)
    ap.add_argument("--min-control-depth", type=float, default=20)
    a = ap.parse_args()
    keep = set(a.species.split(","))
    man = load_manifest(a.manifest)
    bam = pysam.AlignmentFile(a.bam, "rb")

    # per-species control baseline = median depth of single-copy control amplicons
    controls = {}
    for sp in keep:
        d = [region_depth(bam, m["pan_contig"], m["start"], m["end"])
             for m in man if m["role"] == "control" and m["species"] == sp]
        controls[sp] = statistics.median(d) if d else 0.0

    rows = []
    for m in man:
        if m["role"] not in ("cnv", "del") or m["species"] not in keep:
            continue
        sp = m["species"]; base = controls.get(sp, 0.0)
        dep = region_depth(bam, m["pan_contig"], m["start"], m["end"])
        if base < a.min_control_depth:
            call = "no_call(low control depth)"; ratio = ""
        elif m["role"] == "cnv":
            ratio = dep / base
            call = "AMPLIFIED (%.1f copies)" % ratio if ratio >= a.amp_thresh else "single copy"
        else:  # deletion
            ratio = dep / base
            call = "DELETED" if ratio <= a.del_thresh else "present"
        note = CNV_GENES.get(m["gene"]) or DEL_GENES.get(m["gene"]) or ""
        rows.append([a.sample, sp, m["gene"], m["role"], round(dep, 1), round(base, 1),
                     (round(ratio, 2) if ratio != "" else ""), call, note])

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "species", "gene", "type", "depth", "control_baseline", "ratio", "call", "significance"])
        w.writerows(rows)
    flags = [f"{r[2]} {r[7]}" for r in rows if "AMPLIFIED" in r[7] or "DELETED" in r[7]]
    sys.stderr.write(f"[cnv_deletion] {a.sample}: " + ("; ".join(flags) if flags else "no CNV/deletion flags") + "\n")


if __name__ == "__main__":
    main()
