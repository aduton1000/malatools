# References

| File | In git? | Purpose |
|------|---------|---------|
| `pan_reference.fa.gz` | ✅ committed | Combined 6-nuclear + 6-mito, species-tagged reference. **The only genome MalaCall (analysis) needs.** |
| `mito/`, `barcodes/`, `human_*.fasta` | ✅ committed | Small derived/support files. |
| `nuclear/*.fna`, `anno/*.gff` | ❌ git-ignored | Raw per-species genomes + annotations (~250 MB). Only needed to **re-run the MalaForge design**. |

## Getting the raw genomes (design only)

The raw assemblies are large and public, so they are not committed. Download them from NCBI with:

```bash
malaforge fetch-refs            # or: python3 refs/fetch_genomes.py
malaforge fetch-refs --list     # preview the accession table (no network)
malaforge fetch-refs --only Pf3D7,PkHl
```

This writes `refs/nuclear/<name>.fna` and `refs/anno/<name>.gff` for the six taxa.

### Assemblies used

| Name | Accession | Species |
|------|-----------|---------|
| Pf3D7   | GCF_000002765.6 | *P. falciparum* 3D7 |
| PvP01   | GCA_900093555.1 | *P. vivax* P01 |
| PmUG01  | GCF_900090045.1 | *P. malariae* UG01 |
| PkHl    | GCF_000006355.2 | *P. knowlesi* strain H |
| PocGH01 | GCA_900090035.2 | *P. ovale curtisi* GH01 |
| PowCR01 | GCA_900090025.2 | *P. ovale wallikeri* CR01 |

> **Verify:** these accessions reproduce the assemblies used to build the committed panel. Confirm them
> against NCBI before relying on exact coordinates for a fresh design run. Edit the table at the top of
> `fetch_genomes.py` to change them.
