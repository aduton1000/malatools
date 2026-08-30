# Non-falciparum resistance markers — curated (for Module R)

Synthesised from a structured literature review (2026-08-28). **Key finding: the only gene with
functionally validated resistance-associated mutations outside Pf/Pv is `dhfr` (pyrimethamine), and
only in *P. malariae* and *P. ovale*.** Everything else is a *discovery* target (novel/uncharacterised
polymorphism), which the panel amplifies but reports as candidate-grade, not validated.

## Validated / functionally-supported (type + call these)

| Species | Gene | Gene ID | Codons (mutant) | Pf/Pv homolog | Evidence | Note |
|---|---|---|---|---|---|---|
| P. malariae | pmdhfr | PmUG01_05034700 | **R58, N114, F57L** | R58≈PfC59R, N114≈PfS108N | CRISPR-validated (Ibrahim 2024, PMID 39738025); Khim 2012 | **Reference PmUG01 carries the MUTANT allele (R58,N114)** — call vs the SENSITIVE allele S58/S114 (same wt_aa rule as our 7 Pf/Pv markers) |
| P. ovale curtisi | pocdhfr | (PocGH01, tag unverified) | **A15S, S58R, S113N/T** | S58R≈PfC59R, S113N≈PfS108N, F57L≈PfN51I | Epi + E.coli-complementation (Chen 2022; Joste 2024, PMID 38761813) | S58R ~52%, strongest in Central Africa |
| P. ovale wallikeri | powdhfr | (PowCR01, tag unverified) | F57L, S58R, S113N | as above | mostly WT (~10% mutant) | disparate subspecies selection |

## Discovery targets — NO validated marker (amplify for emergence surveillance, report candidate)

- **dhps** (all): Pm/Po/Pk all WT at Pf 436/437/540/581/613; only novel out-of-pocket SNPs.
- **crt** (all): no K76T homolog in Pm/Po/Pk; chloroquine remains clinically effective.
- **mdr1** (all): single-copy; no N86Y/Y184F/D1246Y homolog.
- **k13** propeller (all): WT propeller in Po/Pk; Pm P533L (1 isolate, not validated).
- **cytb** (all): atovaquone Y268-homolog uncharacterised (Pk entirely unstudied — do NOT assume conserved).
- **P. knowlesi**: NO validated marker in any gene; pkdhfr/pkdhps polymorphic but WT at Pf-homologous sites;
  clinically CQ- and ART-sensitive.

## Design implications
- Module R **validated component** = pmdhfr (57/58/114) + podhfr (15/57/58/113) for Poc & Pow, plus the
  full Pf (WHO compendium) + Pv marker set carried over from low-cost-mal-seq.
- Module R **discovery component** = orthologous dhfr/dhps/crt/mdr1/k13-propeller/cytb amplicons across
  all species, reported as candidate-grade with novel-variant scanning.
- Reconcile pmdhfr numbering polarity against the PmUG01 reference before hard-coding (reference = mutant).
