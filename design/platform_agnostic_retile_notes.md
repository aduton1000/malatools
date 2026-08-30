# Platform-agnostic re-tiling — design notes

## Goal
Make the panel **truly platform-agnostic** (one design that runs on ONT, Illumina 2×250, and MGI) and
collapse it toward a **single pool**, with **no compromise on marker coverage or specificity**.

## The problem we found
The original panel was ONT-first: many resistance amplicons were **700–1208 bp** long single reads
covering several codons each. On Illumina 2×250 the middle of a >500 bp amplicon is **unreadable** —
so those markers would be lost. Tightening only the diversity loci would not fix this.

## The rule
**Every pooled amplicon ≤ 400 bp core** (CS1/CS2 tails add 44 → ≤ 444 bp; 2×250 reads then overlap by
≥ ~56 bp for a clean merged consensus). Non-overlapping tiles so they can share one pool.

## What we did
`design/design_retile_resistance.py` re-tiles the **entire** resistance module into ≤400 bp
non-overlapping amplicons, from two sources uniformly:

- **Codon-mapped genes** (all Pf, all Pv, non-Pf dhfr): cluster the validated marker codons
  (`targets/codon_to_genomic_pan.tsv`) into groups each spanning ≤ ~300 bp genomic, **placing every tile
  junction in the largest inter-codon gap** so no marker codon ever lands in a primer/junction blind
  spot. One ≤400 bp amplicon per cluster (e.g. pfcrt's 24 codons → 6 tiles; pfmdr1 → 3; pfdhps → 2).
- **Discovery genes** (non-Pf dhps/k13/crt/mdr1, Pk dhfr): tblastn the Pf protein → map the marker
  window to the species genome → tile sequentially, each tile's template anchored where the previous
  amplicon ends (`design_seq`), which **guarantees strictly non-overlapping** tiles.

Every tile is validated: single product on its own genome, species-specificity by 3′-anchor counts
across the other five taxa, product ≤400 bp.

Other modules were brought to ≤400 too: new-species diversity windows tightened (all 15 now ≤400 **and
all species-specific**); Pk-cluster diagnostics capped (chr1 328 / chr6 341 / chr8 352); three long CNV
amplicons in Module D re-tiled (pfmdr1, pvmdr1, pfpm2); the primary species-ID **cox3** re-designed to
358 bp (see below). Cross-primer dimers were re-engineered by `design/dimer_resolve.py`.

## Result
- **289 pooled amplicons** (S 1 + R 87 + D 15 + V 186), **max 397 bp core / 441 bp tailed — all
  2×250-mergeable.** Resistance grew 43 → 87 (long amplicons became short non-overlapping tiles); **no
  gene and no marker codon was lost.**
- **0 human off-target.**
- **2 pools (286 + 3).** Pool 2 = the pfk13 propeller tiles + one mdr1 tile (see below). The two flagged
  cross-dimers both involve a pool-2 amplicon (a 59.2 °C pair that is *already* separated as overlapping
  k13 tiles — they never share a tube — and one marginal 40.4 °C interaction on a separated amplicon), so
  there are **no functional dimers within either pool.**

## The one irreducible constraint: pfk13
A **true single pool is not achievable** while keeping every amplicon ≤400 bp *and* fully covering the
artemisinin-resistance **k13 propeller**. Its 33 validated/candidate codons are packed into ~786 bp with
essentially **one** usable inter-codon gap in the dense core (476↔493, 51 bp). Covering them all at
≤400 bp therefore forces the two central tiles to **overlap**, and overlapping amplicons must sit in
different pools (co-pooling them lets the outer primers make a chimeric product that steals template and
degrades yield — the exact quality loss we were asked to avoid). So the minimum is:

> **one ~286-plex universal pool + a tiny "k13 propeller" pool (2 tiles).**

This is the *quality-preserving* answer, not a compromise. Dropping k13 candidate codons or using a
>500 bp k13 amplicon (unreadable-in-the-middle on Illumina) were both rejected.

## Species-ID at ≤400 bp
The primary cox3 amplicon was tightened from ~500 bp (12-SNP minimum) to **358 bp (8-SNP minimum**,
Poc↔Pow; Pv↔Pk 10). Eight fixed differences over a 358 bp read is **decisive** for species calling — no
capability loss, and now merge-clean on 2×250. The longer **cox3_long (940 bp, 22 SNP)** and orthogonal
**cox1_S2 (592 bp)** are retained as **OPTIONAL-ONT add-ons** (higher Poc/Pow resolution) that are
**excluded from the universal ≤400 pool** and run separately on ONT when desired.

## Reproduce
`design/design_retile_resistance.py` → `module_R_resistance.tsv` + `retile_overlap_conflicts.tsv`;
`design/design_diversity_newsp.py`, `design/design_pk_cluster.py` (≤400 caps);
`design/dimer_resolve.py`; then `design/assemble_panel.py` (reads the conflicts file and excludes
OPTIONAL-ONT amplicons from the pool).
