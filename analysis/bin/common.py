#!/usr/bin/env python3
"""
common.py — shared helpers for the pan-plasmodium analysis pipeline.
Manifest loading, codon/haplotype pileup from BAM (pysam), translation, depth, IUPAC.
"""
import os, csv
import pysam

BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON2AA = {c1 + c2 + c3: AAS[i * 16 + j * 4 + k]
            for i, c1 in enumerate(BASES) for j, c2 in enumerate(BASES) for k, c3 in enumerate(BASES)}
COMP = str.maketrans("ACGTN", "TGCAN")
MITO_ACC = {"Pf": "NC_002375.1", "Pv": "NC_007243.1", "Poc": "CM190334.1",
            "Pow": "CM190336.1", "Pm": "LT594637.1", "Pk": "NC_007232.1"}


def rc(s):
    return s.translate(COMP)[::-1]


def translate(nt):
    nt = nt.upper()
    return "".join(CODON2AA.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def species_of(gene):
    g = gene.lower()
    if g.startswith("pf"):
        return "Pf"
    if g.startswith("pv"):
        return "Pv"
    if g.startswith("pm"):
        return "Pm"
    if g.startswith("poc"):
        return "Poc"
    if g.startswith("pow"):
        return "Pow"
    if g.startswith("pk"):
        return "Pk"
    return None


def pan_contig(sp, contig):
    """map a species-genome contig name to the pan-reference contig id (sp|contig or sp|MITO)."""
    if contig == MITO_ACC.get(sp):
        return f"{sp}|MITO"
    return f"{sp}|{contig}"


def load_manifest(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows:
        r["start"] = int(r["start"]); r["end"] = int(r["end"]); r["product_bp"] = int(r["product_bp"])
    return rows


def region_depth(bam, contig, start, end, min_mapq=10):
    """mean read depth over [start,end) using pysam."""
    try:
        cov = bam.count_coverage(contig, start, end, quality_threshold=0,
                                 read_callback=lambda r: (not r.is_unmapped) and r.mapping_quality >= min_mapq)
    except (ValueError, KeyError):
        return 0.0
    n = end - start
    if n <= 0:
        return 0.0
    per = [cov[0][i] + cov[1][i] + cov[2][i] + cov[3][i] for i in range(n)]
    return sum(per) / n


def region_read_count(bam, contig, start, end, min_mapq=10):
    try:
        return bam.count(contig, start, end,
                         read_callback=lambda r: (not r.is_unmapped) and r.mapping_quality >= min_mapq)
    except (ValueError, KeyError):
        return 0


def codon_alleles(bam, contig, gstart, gend, strand, min_bq=13, min_mapq=10):
    """Return {aa: fraction, ...} + depth for the codon at 1-based genomic [gstart,gend] (3 bp),
    reading each spanning read's 3 bases, orienting by strand, translating."""
    positions = list(range(gstart - 1, gend))            # 0-based ref positions of the 3 codon bases
    if len(positions) != 3:
        positions = [gstart - 1, gstart, gstart + 1]
    counts = {}
    depth = 0
    try:
        for col_reads in [1]:                            # single pass over reads
            for read in bam.fetch(contig, positions[0], positions[2] + 1):
                if read.is_unmapped or read.mapping_quality < min_mapq or read.is_secondary or read.is_supplementary:
                    continue
                pairs = {ref: q for q, ref in read.get_aligned_pairs(matches_only=True)}  # ref_pos -> read_pos
                rp = read.query_sequence
                qual = read.query_qualities
                bases = []
                ok = True
                for pos in positions:
                    if pos not in pairs:
                        ok = False; break
                    qi = pairs[pos]
                    if qual is not None and qual[qi] < min_bq:
                        ok = False; break
                    bases.append(rp[qi])
                if not ok:
                    continue
                codon = "".join(bases)
                if strand == "-":
                    codon = rc(codon)
                aa = CODON2AA.get(codon.upper(), "X")
                counts[aa] = counts.get(aa, 0) + 1
                depth += 1
    except (ValueError, KeyError):
        return {}, 0
    if depth == 0:
        return {}, 0
    return {aa: c / depth for aa, c in counts.items()}, depth


def haplotype_reads(bam, contig, start, end, min_mapq=10, min_len_frac=0.8, max_reads=4000):
    """Return list of full-length read sequences spanning [start,end) as aligned strings
    (ref-anchored, insertions dropped) for microhaplotype calling."""
    out = []
    span = end - start
    try:
        for read in bam.fetch(contig, start, end):
            if read.is_unmapped or read.mapping_quality < min_mapq or read.is_secondary or read.is_supplementary:
                continue
            pairs = read.get_aligned_pairs(matches_only=True)
            refmap = {rpos: read.query_sequence[qi] for qi, rpos in pairs if rpos is not None}
            covered = [refmap.get(p, "-") for p in range(start, end)]
            if sum(1 for c in covered if c != "-") < span * min_len_frac:
                continue
            out.append("".join(covered))
            if len(out) >= max_reads:
                break
    except (ValueError, KeyError):
        return []
    return out
