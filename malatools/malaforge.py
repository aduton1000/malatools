#!/usr/bin/env python3
"""malaforge — panel-design toolkit for the pan-Plasmodium low-cost amplicon assay.

Designs and validates the species-ID, resistance, diagnostic/CNV and diversity amplicons across all
6 taxa, tails them with the CS1/CS2 linkers, screens cross-dimers, and pools the panel. Pair it with
MalaCall, the analysis pipeline that reads the data this panel generates.
"""
import sys
from ._common import dispatch

TAGLINE = "pan-Plasmodium panel-design toolkit"
EPILOG = ("Most users run the finished panel (see design/panel_tailed.tsv) and order it from the "
          "workbook in order/. These commands reproduce or extend the design. Docs: docs/DESIGN-DOSSIER.")

COMMANDS = [
    ("assemble",       "Tail (CS1/CS2), screen cross-dimers and pool the full panel"),
    ("retile",         "Re-tile resistance markers into ≤400 bp platform-agnostic amplicons"),
    ("species-id",     "Design the cox3 species-ID amplicon on the 6-taxa MSA"),
    ("resistance",     "Design a species-specific resistance amplicon over a marker window"),
    ("multiexon",      "Design an exon-aware amplicon (crt/mdr1) via tblastn codon mapping"),
    ("diversity",      "Design new-species diversity amplicons over conserved-flanked antigens"),
    ("pk-cluster",     "Design the P. knowlesi subpopulation (Mf/Mn/cluster-3) panel"),
    ("orthologs",      "Locate resistance-gene orthologs across the reference set"),
    ("ispcr",          "In-silico PCR / primer-pair validation against the references"),
    ("dimer",          "Re-engineer primers flagged for cross-dimers"),
    ("barcodes",       "Generate the orthogonal 24-nt barcode set"),
    ("check-barcodes", "Check barcodes vs primers, CS tags and the genome for conflicts"),
    ("fetch-refs",     "Download the raw reference genomes for design (NCBI; git-ignored, ~250 MB)"),
]

HANDLERS = {
    "assemble":       "design/assemble_panel.py",
    "retile":         "design/design_retile_resistance.py",
    "species-id":     "design/design_species_id.py",
    "resistance":     "design/design_resistance_amplicon.py",
    "multiexon":      "design/design_multiexon_amplicon.py",
    "diversity":      "design/design_diversity_newsp.py",
    "pk-cluster":     "design/design_pk_cluster.py",
    "orthologs":      "design/locate_orthologs.py",
    "ispcr":          "design/in_silico_pcr.py",
    "dimer":          "design/dimer_resolve.py",
    "barcodes":       "analysis/bin/generate_barcodes.py",
    "check-barcodes": "analysis/bin/check_barcode_conflicts.py",
    "fetch-refs":     "refs/fetch_genomes.py",
}


def main():
    sys.exit(dispatch("malaforge", TAGLINE, COMMANDS, HANDLERS, EPILOG))


if __name__ == "__main__":
    main()
