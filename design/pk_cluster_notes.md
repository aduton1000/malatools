# P. knowlesi cluster-typing add-on — design notes

## Purpose
*P. knowlesi* is not one population. It splits into three deeply divergent subpopulations that differ in
macaque host, geography, and (potentially) clinical behaviour:

| cluster | name | macaque host | geography |
|---|---|---|---|
| 1 | **Mf** | *Macaca fascicularis* (long-tailed) | Malaysian Borneo (~⅔ of Borneo human cases) |
| 2 | **Mn** | *Macaca nemestrina* (pig-tailed) | Malaysian Borneo (~⅓) |
| 3 | — | (mainland) | Peninsular Malaysia; includes the lab clone H / reference |

The Mf and Mn nuclear genomes are dimorphic, with a large **Mn→Mf introgression on chromosome 8**
(mosquito/sexual-stage genes). None of this is visible to the species-ID or antigen-diversity modules —
this add-on adds it.

## Two tiers (8 amplicons), all primary-source-grounded

Everything is anchored by **primer/gene sequence, not by published coordinates.** The Divis papers report
positions in the PKNH (GeneDB) assembly, which is **offset from the RefSeq strain-H assembly
(GCF_000006355.2)** used across this project — e.g. the NC03_2 locus sits ~45 kb from its printed
coordinate, and one chr12 locus >1 Mb away. Sequence anchoring makes the design assembly-version-proof.

### Tier 1 — within/between-cluster diversity (5 microsatellites)
The microsatellite panel of **Divis et al. 2015** (PLoS Pathog, PMID 25993697, S1 Table) is the de-facto
*P. knowlesi* diversity standard; genotyped across isolates it defines Mf vs Mn and gives
relatedness/COI. We take each locus' **nest-2 primer pair verbatim** and keep those that still prime the
clone-H reference. **6 of the 10 locate cleanly** (≤3 mismatch, correct product size); one (CD11_157)
falls just outside our exact-3′-seed rule, leaving **5 carried**:

| amplicon | chr | product | note |
|---|---|---|---|
| pkclu_NC03_2  | 3  | 154 bp | **FST > 0.3 — top Mf/Mn differentiator** |
| pkclu_CD08_61 | 8  | 231 bp | inside the chr8 **cap380** introgression window |
| pkclu_NC10_1  | 10 | 275 bp | |
| pkclu_NC12_2  | 12 | 347 bp | |
| pkclu_NC12_4  | 12 | 235 bp | |

The other 4 (CD05_06, NC09_1, CD13_61, CD13_107) have **diverged in clone H** (5–8 primer mismatches) —
their capillary primers were built on Borneo isolates and cannot be validated against this reference; they
are dropped honestly rather than forced.

### Tier 2 — direct cluster-diagnostic SNP amplicons (3)
These are the **cluster callers**. Rebuilt as **sequencing-native** amplicons: instead of the published
allele-specific PCR (band present/absent), we place primers that amplify the locus in every sample and put
the fixed diagnostic SNP **inside the read**, so the sequencer calls the allele directly.

| amplicon | chr | product | axis | basis |
|---|---|---|---|---|
| pkclu_diag_chr1_C1A    | 1 | 385 bp | **Mf vs Mn** | Divis 2020 C1A/C2F fixed-difference locus (EID 10.3201/eid2607.190924); anchored on shared reverse CX01AR (0-mismatch on clone H) |
| pkclu_diag_chr6_C2J    | 6 | 544 bp | **Mn**       | Divis 2020 C2J fixed-difference locus; anchored on C206DR (0-mismatch) |
| pkclu_diag_chr8_cap380 | 8 | 352 bp | **introgression** | cap380 = PKNH_0820800 (933.1–943.3 kb), the introgression's strongest gene (~101 SNPs FST>0.4; Assefa 2015 PMID 26438871 / Diez Benavente 2017 10.1371/journal.pgen.1007008); anchored on the RefSeq gene annotation |

The C1A/C2J loci come from the **9,293 SNPs with complete fixed difference (FST = 1)** between Mf and Mn.
The chr8 cap380 amplicon reads the Mn-type-vs-Mf-type allele that marks the introgression (~73 % of
Betong Mf isolates carry the introgressed Mn cap380 allele).

### Cluster-3 (Peninsular) axis
Read from the **organellar** signal already in the panel: the mito cox3/cox1 ID amplicons' haplotypes
separate the Borneo vs Peninsular lineages (Diez Benavente 2017: organellar SNPs resolve three haplotype
groups; the Peninsular/cultured clade is distinct). No extra amplicon needed — a documented reuse.

## How a sample is typed
1. **cap380 + C1A + C6/C2J** read the Mf-vs-Mn nuclear axis and the chr8 introgression allele → Mf / Mn / mixed.
2. **Organellar (cox3/cox1)** haplotype → Borneo vs **Peninsular (cluster 3)**.
3. **5 microsatellites** → within-cluster relatedness / COI (allele = repeat length + flanking SNP haplotype, read directly from sequence).

## Honest findings & limitations
1. **P. knowlesi is intrinsically microsatellite-poor** — 0 perfect trinucleotide runs ≥6 copies
   genome-wide (vs hundreds in *P. falciparum*). This is why the Divis alleles are short in clone H and why
   de-novo STR discovery is unproductive; the diagnostic SNP amplicons (Tier 2) therefore carry the primary
   cluster signal, with the 5 microsatellites for relatedness.
2. **Only 5/10 Divis loci re-prime the clone-H reference.** The full 10-locus panel would need the PKNH
   assembly (or per-locus primer redesign on Borneo-isolate sequence) — a wet-lab/assembly task.
3. **Reference is cluster 3.** The clone-H genome is Peninsular, so it is one allele at every locus; the
   fixed-difference calls are validated against the published SNP sets, not re-derived here from a
   population VCF. A denser barcode could be mined from the ~1.6 M-SNP WGS catalogues (Assefa 2015;
   Diez Benavente 2017) — future work.
4. **Coordinate provenance:** published positions are PKNH; we anchor by sequence and report clone-H
   (GCF_000006355.2) positions. Confirm concordance before any coordinate-based downstream step.

## Reproduce
`design/design_pk_cluster.py` → `design/module_V_pk_cluster.tsv`; folded into Module V by
`design/assemble_panel.py`. Divis primer sequences transcribed from PMID 25993697 S1 Table and
10.3201/eid2607.190924 Appendix Table 1; cap380 = PKNH_0820800 from the RefSeq annotation `PkHl.gff`.
