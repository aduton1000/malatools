#!/usr/bin/env python3
"""
simulate.py — generate known-truth synthetic reads for the whole panel and the matching truth tables,
so the pipeline can be validated end-to-end against ground truth.

Reads are built directly from the pan-reference at each amplicon's mapped coordinates (forward strand),
then edited to encode the truth: resistance codons, CNV (read-count scaling), gene deletions (drop reads),
COI (k microhaplotypes per diversity locus), Pk cluster (diagnostic-SNP allele), plus a substitution error
model. Single-end full-amplicon reads (ONT-style) — the panel's native mode; also the strictest test of
mapping + pileup.

Usage: simulate.py <manifest.tsv> <pan_reference.fa.gz> <codon_map.tsv> <pk_panel.tsv> <outdir>
Writes: <outdir>/reads/<sample>.fastq.gz and <outdir>/truth/*.tsv
"""
import sys, os, csv, gzip, argparse, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../bin"))
from common import CODON2AA, rc, species_of, pan_contig, load_manifest

AA2CODON = {}
for cod, aa in CODON2AA.items():
    AA2CODON.setdefault(aa, cod)
BASE_DEPTH = 220
ERR = 0.006
random.seed(42)

# ---- truth definition: a comprehensive panel exercising every output type ----
SAMPLES = {
    "S1_Pf_multidrug": {
        "species": {"Pf": 1.0},
        "resistance": {"Pf": {"pfcrt": {76: "T"}, "pfdhfr": {108: "N", 51: "I", 59: "R"},
                              "pfdhps": {437: "G", 540: "E"}, "pfk13": {580: "Y"}, "pfmdr1": {86: "Y"}}},
        "cnv": {"Pf": {"pfmdr1": 3}},
        "deletion": {"Pf": {"pfhrp2": True, "pfhrp3": False}},
        "coi": {"Pf": 2},
    },
    "S2_Pf_sensitive": {
        "species": {"Pf": 1.0}, "resistance": {}, "cnv": {}, "deletion": {}, "coi": {"Pf": 1},
    },
    "S3_PfPv_mixed": {
        "species": {"Pf": 0.65, "Pv": 0.35},
        "resistance": {"Pf": {"pfcrt": {76: "T"}}, "Pv": {"pvdhfr": {58: "R", 117: "N"}}},
        "cnv": {}, "deletion": {}, "coi": {"Pf": 1, "Pv": 2},
    },
    "S4_Pm_resistant": {   # PmUG01 reference CARRIES the mutant at 58/114 -> tests polarity (truth = mutant)
        "species": {"Pm": 1.0}, "resistance": {"Pm": {"pmdhfr": {58: "R", 114: "N"}}}, "coi": {"Pm": 1},
    },
    "S5_Pm_sensitive": {   # truth = sensitive S/S even though the reference carries R/N -> hardest polarity test
        "species": {"Pm": 1.0}, "resistance": {"Pm": {"pmdhfr": {58: "S", 114: "S"}}}, "coi": {"Pm": 1},
    },
    "S6_Pk_Mn": {"species": {"Pk": 1.0}, "resistance": {}, "coi": {"Pk": 1}, "pk_cluster": "Mn"},
    "S7_PocPow_mixed": {"species": {"Poc": 0.5, "Pow": 0.5}, "resistance": {}, "coi": {"Poc": 1, "Pow": 1}},
    "S8_negative": {"species": {}, "resistance": {}, "coi": {}},   # no parasite
    # a clonally-related pair (same hap_seed -> shared lineages) to validate relatedness; both Pf monoclonal
    "S9_Pf_clonalA": {"species": {"Pf": 1.0}, "resistance": {}, "coi": {"Pf": 1}, "hap_seed": "CLONE_K"},
    "S10_Pf_clonalB": {"species": {"Pf": 1.0}, "resistance": {}, "coi": {"Pf": 1}, "hap_seed": "CLONE_K"},
    # extra unrelated Pv samples so Pv has a population (COI-MLE / allele freqs need >1 sample per species)
    "S11_Pv": {"species": {"Pv": 1.0}, "resistance": {}, "coi": {"Pv": 1}},
    "S12_Pv": {"species": {"Pv": 1.0}, "resistance": {}, "coi": {"Pv": 1}},
    "S13_Pv": {"species": {"Pv": 1.0}, "resistance": {}, "coi": {"Pv": 2}},
}
RELATED_TRUTH = [("S9_Pf_clonalA", "S10_Pf_clonalB")]   # expected highly related


def encode_codon(gene_aa, strand):
    cod = AA2CODON.get(gene_aa, "NNN")
    return rc(cod) if strand == "-" else cod        # write onto the forward strand


def load_needed_refs(refgz, contigs):
    want = set(contigs); ref = {}; name = None; buf = []
    with gzip.open(refgz, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name in want:
                    ref[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            elif name in want:
                buf.append(line.strip())
        if name in want:
            ref[name] = "".join(buf)
    return ref


def load_codons(path):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        try:
            r["codon"] = int(r["codon"]); r["genomic_start"] = int(r["genomic_start"])
        except (ValueError, TypeError):
            continue
        out.append(r)
    return out


def mutate_read(seq, err):
    if err <= 0:
        return seq
    b = list(seq); alt = "ACGT"
    for i in range(len(b)):
        if random.random() < err:
            b[i] = random.choice([x for x in alt if x != b[i].upper()])
    return "".join(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest"); ap.add_argument("ref"); ap.add_argument("codon_map")
    ap.add_argument("pk_panel"); ap.add_argument("outdir")
    a = ap.parse_args()
    man = load_manifest(a.manifest)
    refs = load_needed_refs(a.ref, {m["pan_contig"] for m in man})
    codons = load_codons(a.codon_map)
    # index codons by (species, gene); add a bare-gene alias so manifest genes like "dhfr" (species Pm)
    # match codon-map genes like "pmdhfr".
    cmap = {}
    PFX = {"Pf": "pf", "Pv": "pv", "Pm": "pm", "Poc": "poc", "Pow": "pow", "Pk": "pk"}
    for c in codons:
        sp = species_of(c["gene"]); g = c["gene"]
        cmap.setdefault((sp, g), []).append(c)
        pfx = PFX.get(sp, "")
        if pfx and g.startswith(pfx):
            cmap.setdefault((sp, g[len(pfx):]), []).append(c)   # alias: pmdhfr -> ("Pm","dhfr")
    pk_panel = {r["amplicon"]: r for r in csv.DictReader(open(a.pk_panel), delimiter="\t")} if os.path.exists(a.pk_panel) else {}

    os.makedirs(f"{a.outdir}/reads", exist_ok=True); os.makedirs(f"{a.outdir}/truth", exist_ok=True)
    truth_rows = {"species": [], "resistance": [], "cnv": [], "deletion": [], "coi": [], "cluster": [],
                  "related": [list(p) for p in RELATED_TRUTH]}

    for sname, spec in SAMPLES.items():
        sp_frac = spec.get("species", {})
        for sp in ["Pf", "Pv", "Poc", "Pow", "Pm", "Pk"]:
            truth_rows["species"].append([sname, sp, "present" if sp_frac.get(sp, 0) > 0 else "absent"])
        for sp, genes in spec.get("resistance", {}).items():
            for gene, cod in genes.items():
                for codon, aa in cod.items():
                    truth_rows["resistance"].append([sname, sp, gene, codon, aa])
        for sp, g in spec.get("cnv", {}).items():
            for gene, cn in g.items():
                truth_rows["cnv"].append([sname, sp, gene, cn])
        for sp, g in spec.get("deletion", {}).items():
            for gene, deld in g.items():
                truth_rows["deletion"].append([sname, sp, gene, "DELETED" if deld else "present"])
        for sp, k in spec.get("coi", {}).items():
            truth_rows["coi"].append([sname, sp, k])
        if "pk_cluster" in spec:
            truth_rows["cluster"].append([sname, spec["pk_cluster"]])

        reads = []
        rid = 0
        for m in man:
            sp = m["species"]
            frac = sp_frac.get(sp, 0.0)
            if frac <= 0:
                continue
            template = refs[m["pan_contig"]][m["start"]:m["end"]]
            if not template:
                continue
            gene = m["gene"]
            # deletion: drop reads for a deleted gene
            deld = spec.get("deletion", {}).get(sp, {}).get(gene, False)
            if deld:
                continue
            # base read count, scaled by species fraction and CNV copy number
            copies = spec.get("cnv", {}).get(sp, {}).get(gene, 1)
            ndepth = max(1, int(BASE_DEPTH * frac * copies))

            # build haplotype templates
            haps = [list(template)]
            # resistance edits (all clones carry the truth allele)
            for c in cmap.get((sp, gene), []):
                # truth is keyed by the codon-map gene name (e.g. "pmdhfr"), not the manifest gene ("dhfr")
                aa = spec.get("resistance", {}).get(sp, {}).get(c["gene"], {}).get(c["codon"])
                if aa is None:
                    continue
                off = (c["genomic_start"] - 1) - m["start"]
                if 0 <= off <= len(template) - 3:
                    codon = encode_codon(aa, c["strand"])
                    for h in haps:
                        h[off], h[off + 1], h[off + 2] = codon[0], codon[1], codon[2]
            # Pk cluster diagnostic allele
            if m["amplicon"] in pk_panel and "pk_cluster" in spec:
                pr = pk_panel[m["amplicon"]]; off = int(pr["offset"]); allele = pr.get(spec["pk_cluster"])
                if allele and 0 <= off < len(template):
                    for h in haps:
                        h[off] = allele
            # COI: diversity loci get k SAMPLE-SPECIFIC microhaplotypes (so unrelated samples differ and
            # a related pair sharing hap_seed shares lineages -> lets relatedness be validated).
            k = spec.get("coi", {}).get(sp, 1)
            if m["role"] == "diversity" and len(template) > 40:
                import zlib
                seed = spec.get("hap_seed", sname)
                base = list(template)
                nvar = min(5, (len(base) - 20) // 3)
                lr = random.Random(zlib.crc32(f"{seed}|{m['amplicon']}".encode()))
                vpos = lr.sample(range(10, len(base) - 10), nvar)
                haps = []
                for ci in range(k):
                    hh = list(base)
                    cr = random.Random(zlib.crc32(f"{seed}|{m['amplicon']}|{ci}".encode()))
                    for vp in vpos:
                        hh[vp] = cr.choice("ACGT")
                    haps.append(hh)
            # emit reads across haplotypes
            for i in range(ndepth):
                h = haps[i % len(haps)]
                seq = mutate_read("".join(h).upper(), ERR)
                reads.append((f"{sname}_{m['amplicon']}_{rid}", seq)); rid += 1

        with gzip.open(f"{a.outdir}/reads/{sname}.fastq.gz", "wt") as fh:
            for name, seq in reads:
                fh.write(f"@{name}\n{seq}\n+\n{'I'*len(seq)}\n")
        sys.stderr.write(f"[simulate] {sname}: {len(reads)} reads\n")

    for key, rows in truth_rows.items():
        with open(f"{a.outdir}/truth/{key}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t"); w.writerows(rows)
    print(f"[simulate] wrote {len(SAMPLES)} samples + truth tables -> {a.outdir}")


if __name__ == "__main__":
    main()
