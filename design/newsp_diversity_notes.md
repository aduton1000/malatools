# New-species diversity panel — design notes (Pm / Poc / Pow / Pk)

## Why this is different from the Pf/Pv diversity panel

The Pf and Pv microhaplotype panels in Module V (163 amplicons) were selected from **population
genomics** — thousands of sequenced isolates identify neutral, SNP-dense, high-heterozygosity loci that
resolve complexity-of-infection (COI), relatedness, and structure. *P. malariae*, *P. ovale curtisi*,
*P. ovale wallikeri* and *P. knowlesi* have **no comparable resource**: there is no published SNP
barcode or microhaplotype panel for any of the four, and genome-wide nucleotide diversity is 1–2 orders
of magnitude below sympatric Pf (Pm π ≈ 1.7×10⁻⁴; ovale π ≈ 10⁻⁴; several Pk antigens are under
*purifying* selection). Most historically "diverse" loci in these species are **tandem-repeat / indel**
markers — excellent for capillary size-typing, poor for clean microhaplotype/SNP calling.

We therefore build a **first-pass, homology-anchored diversity panel**: amplicons over the
conserved-flanked, SNP-dense domains of **single-copy antigen genes** that are documented as polymorphic
(balancing / diversifying selection) in the target species or their close relatives. Loci are placed by
tblastn HSP codon→genomic mapping (the same paralog-safe geometry used for resistance), primer3-designed
in the conserved flanks, and validated by in-silico PCR for species-specificity and human off-target.

**These are explicitly discovery-grade.** Locus choice is by homology to known-diverse domains, *not* by
measured intra-species polymorphism (which the single reference genome cannot show). A production panel
would add dozens of neutral windows mined from the published WGS SNP catalogues (below).

## Loci selected (15 amplicons: Pm 4, Poc 4, Pow 4, Pk 3)

| Locus | Pf ortholog | Variable domain (window) | Taxa | Evidence |
|---|---|---|---|---|
| **ama1** | PF3D7_1133400 | Domain I polymorphic cluster (aa 196–300) | Pm, Poc, Pow | Pm highest-π diversifying antigen (PMID 29321885); balancing in Pf/Pv |
| **trap** | PF3D7_1335900 | vWA adhesive core (aa 60–200) | Pm, Pk, Poc, Pow | Pm diversifying, 17 aa variants (PMID 29321885); Pk 74 SNPs / 29 haplotypes (10.1186/s12936-018-2423-1) |
| **csp** | PF3D7_0304600 | C-terminal Th2R/Th3R (aa 300–375) | Pm, Pk, Poc, Pow | Pm 6/9 haplotypes, Fst 0.40 (PMID 30012172); ovale dimorphic + Th2R/Th3R (PMID 31170213); Pk π 0.029 (PMID 26379157) |
| **msp1** | PF3D7_0930300 | C-terminal p19/p42 double-EGF (aa 1620–1710) | Pm, Pk, Poc, Pow | Pk MSP1-42 18 haplotypes, balancing in Borneo (PMID 30518419); Pm block-IV clean SNP window nearby (PMID 36114242) |

### Deliberate exclusions
- **PkAMA1 excluded** — *P. knowlesi* AMA1 Domain I is under **purifying** selection (low diversity,
  26–27 haplotypes; 10.1371/journal.pone.0124400). A weak marker; not designed.
- **CelTOS dropped** — no measured diversity in any of the four target species; it was a pure homology
  guess and adds amplicon load without evidenced payoff.
- **EBA / DBP / NBP / MSP7 families avoided** — lineage-expanded multigene families with **no clean 1:1
  Pf ortholog**; tblastn would anchor ambiguously across paralogs. Excluded to keep locus placement
  paralog-safe. (This forgoes PkDBPαII — a genuinely informative Pk locus — as a deliberate,
  documented trade-off; see future work.)
- **Repeat cores avoided** — CSP central tetrapeptide repeat, MSP1 block-2 / block-8 VNTR, TRAP TSP1
  12-bp repeat, potra/ctrp repeats: all excluded in favour of the flanking SNP-dense windows.

## Specificity
13 / 15 amplicons are species-specific by in-silico PCR. 2 are shared between the ovale sisters
(**Poc↔Pow**: trap, msp1) — the same single case the resistance panel hits, because Poc and Pow are only
~12 mito-SNP apart. Attribution: the cox3 species-ID amplicon resolves Poc vs Pow, so a mono-ovale
infection is correctly assigned; only a *Poc + Pow co-infection* cannot self-attribute these two shared
diversity amplicons (a documented, narrow limitation).

## Honest limitations (carry into the dossier)
1. **Discovery-grade, not validated.** Diversity of each amplicon is inferred from homology to
   known-diverse domains, not measured. Expected segregating SNPs per amplicon are **few** (low π
   everywhere) — COI/relatedness power per locus is lower than a Pf microhaplotype.
2. **Panel is deliberately small (15).** Adequate power for these low-diversity species needs *many*
   loci. This is a seed set; scale-up mines the WGS SNP catalogues:
   - *P. malariae* — 131,601 core-genome SNPs, 251 isolates (10.1038/s41467-024-55102-3).
   - *P. ovale* (both) — 73,015 Poc / 45,669 Pow SNPs (10.1038/s41467-024-54667-3).
   - *P. knowlesi* — ~1.6 M SNPs (PMID 26438871; 10.1371/journal.pgen.1007008).
3. **P. ovale curtisi** has essentially **no prior diversity resource**; wallikeri has only an
   8-microsatellite panel (PMID 37329228). Both are near-universally monoclonal in surveys → limited COI
   signal expected.
4. **P. knowlesi subpopulation structure** (Mf / Mn / cluster-3 macaque-host-associated lineages, defined
   by a chromosome-8 introgression and ~9,293 fixed SNPs; PMID 26438871, 10.3201/eid2607.190924) is **not**
   captured by these antigen amplicons. A dedicated Pk cluster-typing add-on (chr8 fixed-SNP window +
   the Divis-2015 10-microsatellite regions, PMID 25993697) is future work.
5. Antigen diversity data for Pm/ovale are **Asia-biased**; African-population behaviour of these windows
   is unconfirmed.

## Reproduce
`design/design_diversity_newsp.py` (antigen list + windows + per-taxon exclusions at the top) →
`design/module_V_newsp_diversity.tsv` → folded into Module V by `design/assemble_panel.py`.
