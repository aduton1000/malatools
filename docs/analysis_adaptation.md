# Pan-species analysis — how the pipeline extends low-cost-mal-seq

The pan-Plasmodium panel reuses the validated low-cost-mal-seq Nextflow pipeline (demux → align →
resistance → CNV/hrp → novel-variant → interpret → report → cohort), with three additions for
multi-species calling. All per-sample steps stay species-agnostic because they are driven by the
combined reference and the pan codon map.

## 1. Combined reference
`refs/pan_reference.fa` = the six nuclear genomes (Pf3D7, PvP01, PmUG01, PkH, PocGH01, PowCR01) plus the
six mitochondrial genomes, concatenated with species-tagged contig names. minimap2 aligns every read to
this single reference.

## 2. Species detection & identification (new)
- The **cox3 species-ID amplicon** (Module S) reads out species by sequence. In practice this is a
  read tally: reads from the cox3 region map to their species' mitochondrion, so the **per-species cox3
  read count = the species call**, and multiple species above a minor-fraction threshold = a **mixed
  infection**, with read proportions giving approximate per-species density.
- `targets/cox3_species_reference.fasta` holds the six reference inserts for a direct, alignment-free
  confirmation (nearest-reference by SNP distance; ≥12-SNP separation between all taxa, incl. Poc/Pow
  and Pv/Pk).

## 3. Per-species resistance (species-resolved by construction)
- Resistance amplicons are **species-specific** (Module R), so in a mixed infection each species'
  resistance reads come from its own dedicated amplicon and map to that species' reference — resistance
  is attributed to the correct species automatically.
- `resistance_genotyper.py` runs against `targets/codon_to_genomic_pan.tsv` (126 codons: Pf WHO
  compendium + Pv + validated new-species dhfr). It calls each codon **against the sensitive allele
  (`wt_aa`), never the reference** — essential because **9 markers carry the mutant allele in their
  reference genome**, including **pmdhfr 58 (R) and 114 (N)** (PmUG01 is antifolate-mutant), all flagged
  `REF_CARRIES_NONWT`.
- The compendium interpretation, novel-variant discovery, CNV/hrp2-3, diversity (COI/relatedness/
  structure), per-patient reports, and cohort dashboard all reuse the low-cost-mal-seq scripts unchanged.

## 4. Honest reporting
- Non-Pf resistance beyond **dhfr** (dhps/crt/mdr1/k13/cytb) is **discovery-grade** (no validated
  markers) — reported as candidate, with novel-variant scanning for emergence.
- Species proportions are approximate (PCR-efficiency bias); minor-fraction species are flagged, not
  quantified.
- Poc vs Pow are separated by cox3; a rare Poc+Pow co-infection is the only case where the
  ovale-shared resistance amplicons cannot self-attribute.
