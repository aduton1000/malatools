#!/usr/bin/env python3
"""
build_order_workbook.py — generate a vendor-ready procurement workbook for a full pan-Plasmodium run.
Mirrors the Pf/Pv PILOT-ORDER-WORKBOOK structure, extended for 6 taxa / 289-amplicon single-pool panel.

Run size assumption: ONE full run of 96 samples (one 96-well plate / one flow cell). Oligos (R1 primers +
R2 combinatorial barcodes) are ordered ONCE and reused across runs; reagents/flow cells/controls are per run.
"""
import csv, os
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CS1 = "ACACTGACGACATGGTTCTACA"          # Fluidigm CS1 universal linker (external)
CS2 = "TACGGTAGCAGAGACTTGGTCT"          # Fluidigm CS2 universal linker (external)
P5  = "AATGATACGGCGACCACCGAGATCTACAC"   # Illumina P5 (29 nt)
P7  = "CAAGCAGAAGACGGCATACGAGAT"        # Illumina P7 (24 nt)

# ---- load design ----
panel = list(csv.DictReader(open(f"{ROOT}/design/panel_tailed.tsv"), delimiter="\t"))
pools = {r["amplicon"]: r["pool"] for r in csv.DictReader(open(f"{ROOT}/design/panel_pools.tsv"), delimiter="\t")}
bcs   = [(r["barcode_id"], r["sequence"]) for r in csv.DictReader(open(f"{ROOT}/refs/barcodes/barcodes_96.tsv"), delimiter="\t")]
bcmap = dict(bcs)

MOD_LABEL = {"S": "species-ID", "R": "resistance", "D": "diagnostic/CNV", "V": "diversity"}

# ---- styling ----
INK   = "231D19"; PAPER = "FBF8F2"; BLOOD = "9E2B25"; TEAL = "0E6E66"; BAND = "EFEAE0"; GOLD = "8F6B1C"
f_title  = Font(name="Calibri", size=13, bold=True, color=BLOOD)
f_sub    = Font(name="Calibri", size=10, italic=True, color="6C625A")
f_hdr    = Font(name="Calibri", size=10, bold=True, color=PAPER)
f_cell   = Font(name="Calibri", size=10, color=INK)
f_mono   = Font(name="Consolas", size=9, color=INK)
f_note   = Font(name="Calibri", size=9, italic=True, color=GOLD)
fill_hdr = PatternFill("solid", fgColor=INK)
fill_band= PatternFill("solid", fgColor=BAND)
fill_p2  = PatternFill("solid", fgColor="F6E7C9")     # pool-2 highlight
thin     = Side(style="thin", color="DED4C5")
border   = Border(bottom=thin)
wrap     = Alignment(vertical="top", wrap_text=True)
top      = Alignment(vertical="top")


def style_header(ws, row, ncol):
    for c in range(1, ncol + 1):
        cell = ws.cell(row=row, column=c)
        cell.font = f_hdr; cell.fill = fill_hdr; cell.alignment = top


def sheet(wb, name, title, subtitle, headers, widths):
    ws = wb.create_sheet(name)
    ws.cell(1, 1, title).font = f_title
    ws.cell(2, 1, subtitle).font = f_sub
    hr = 4
    for j, h in enumerate(headers, 1):
        ws.cell(hr, j, h)
    style_header(ws, hr, len(headers))
    ws.freeze_panes = ws.cell(hr + 1, 1)
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    return ws, hr + 1


def wells_colmajor():
    return [f"{r}{c}" for c in range(1, 13) for r in "ABCDEFGH"]


# ======================================================================= READ ME
wb = openpyxl.Workbook()
rm = wb.active; rm.title = "READ ME"
rm.column_dimensions["A"].width = 118
lines = [
 ("MalaForge — PANEL ORDER WORKBOOK", f_title),
 ("One low-cost amplicon assay for all 5 human malaria parasites (6 taxa: Pf, Pv, P. ovale curtisi, "
  "P. ovale wallikeri, Pm, Pk). Design: MalaForge · Analysis: MalaCall.", f_sub),
 ("", f_cell),
 ("SCOPE — everything needed to run ONE full 96-sample run of the complete 289-amplicon panel, on any of the "
  "three platforms (ONT / Illumina / MGI), plus controls.", f_cell),
 ("", f_cell),
 ("RUN SIZE: 96 samples (one 96-well plate / one flow cell). To scale: the Round-1 primers and the R2 barcode "
  "oligos are a ONE-TIME purchase and cover every future run; only reagents, flow cells and controls scale "
  "per run. For >96 samples add barcode pairs (the 96-set supports up to 96×96 = 9,216 combinatorially).", f_cell),
 ("", f_cell),
 ("HOW TO USE:", Font(bold=True, size=10, color=INK)),
 (" • Oligo sheets are vendor-ready (Plate / Well / Name / Sequence / Scale / Purification) — hand to any oligo "
  "vendor (IDT, Eurofins, Sigma-Merck, Macrogen).", f_cell),
 (" • Round-1 primers (Oligos-R1-*) are ORDERED ONCE and reused across all runs and all platforms — the panel "
  "is platform-agnostic, so there is a SINGLE universal Round-1 primer set (no separate ONT vs short-read "
  "versions, unlike the Pf/Pv pilot).", f_cell),
 (" • Round-2 barcoding uses COMBINATORIAL dual indexing: 12 forward + 8 reverse barcode oligos give 96 unique "
  "sample tags from just 20 oligos per platform. See 'Barcode-layout-96' for the sample→(fwd,rev) map "
  "(this file is also the demux --barcode_map input).", f_cell),
 (" • Reagent sheets list kits/consumables per platform with indicative catalog numbers + vendor URLs and "
  "quantities sized for 96 samples.", f_cell),
 ("", f_cell),
 ("IMPORTANT — VERIFY: catalog numbers change. Every Cat# is INDICATIVE and must be confirmed against the "
  "vendor's current catalogue before purchase (the 'Verify' column flags this). Oligo/barcode SEQUENCES are "
  "exact.", Font(italic=True, size=10, color=BLOOD)),
 ("", f_cell),
 ("PANEL AT A GLANCE — 289 amplicons, all ≤400 bp core (platform-agnostic; ONT + Illumina 2×250 + MGI):",
  Font(bold=True, size=10, color=TEAL)),
 ("   • Module S (species-ID, cox3)        1 amplicon", f_mono),
 ("   • Module R (resistance)             87 amplicons   (Pf 26, Pk 15, Pm 12, Poc 12, Pow 12, Pv 7, + shared)", f_mono),
 ("   • Module D (diagnostic / CNV)       15 amplicons   (pfhrp2/3 deletion, pfmdr1/pm2-3/gch1/pvmdr1 CNV + controls)", f_mono),
 ("   • Module V (diversity / COI / IBD) 186 amplicons   (Pf 80, Pv 83, new-species 15, Pk cluster 8)", f_mono),
 ("   → 578 Round-1 tailed primers (one-time). 2 pools: Pool 1 = 286 amplicons; Pool 2 = 3 (irreducible pfk13 "
  "propeller pair + one marginal-dimer Pow_mdr1 tile). Pool 2 is a tiny separate 3-plex reaction per sample.",
  f_cell),
 ("", f_cell),
 ("TAILS & BARCODES: Round-1 primers carry Fluidigm CS1/CS2 universal linkers (CS1=" + CS1 + ", CS2=" + CS2 +
  "). Barcodes are our own orthogonal 24-nt set (min Hamming 10) — NOT ONT's licensed kit. Demux with "
  "`malacall demux` (CS1/CS2-anchored).", f_cell),
 ("", f_cell),
 ("CORE run (species+resistance+CNV only) = Oligos-R1-core + one R2 barcode sheet + Reagents-Common + one "
  "platform reagent sheet + Controls. Add Oligos-R1-diversity for the full COI/relatedness panel.", f_cell),
]
for i, (txt, fnt) in enumerate(lines, 1):
    c = rm.cell(i, 1, txt); c.font = fnt; c.alignment = Alignment(wrap_text=True, vertical="top")

# ======================================================================= R1 OLIGO SHEETS
R1_HDR = ["Plate", "Well", "Oligo_Name", "Sequence_5to3", "Length_nt", "Scale", "Purification", "Pool", "Module", "Notes"]
R1_W   = [16, 7, 30, 52, 10, 10, 20, 6, 14, 30]


def emit_r1(ws, start_row, amplicons, plate_prefix):
    wells = wells_colmajor(); wi = 0; plate = 1; r = start_row; n = 0
    for amp in amplicons:
        name = amp["amplicon"]; pool = pools[name]; mod = MOD_LABEL[amp["module"]]
        for suf, seqcol in (("F", "tailed_F(CS1+F)"), ("R", "tailed_R(CS2+R)")):
            if wi >= 96:
                wi = 0; plate += 1
            seq = amp[seqcol]
            note = ""
            if pool == "2":
                note = "POOL 2 — separate 3-plex reaction"
            if any(ch not in "ACGT" for ch in amp["core_F"] + amp["core_R"]):
                note = (note + "; " if note else "") + "degenerate (IUPAC) — hand-mix"
            row = [f"{plate_prefix}_p{plate}", wells[wi], f"{name}.{suf}", seq, len(seq),
                   "25 nmol", "Standard desalting", pool, mod, note]
            for j, v in enumerate(row, 1):
                cell = ws.cell(r, j, v)
                cell.font = f_mono if j in (2, 3, 4) else f_cell
                cell.alignment = wrap if j == 10 else top
                cell.border = border
                if pool == "2":
                    cell.fill = fill_p2
            wi += 1; r += 1; n += 1
    return r, n


core = [a for a in panel if a["module"] in ("S", "R", "D")]
div  = [a for a in panel if a["module"] == "V"]

ws, sr = sheet(wb, "Oligos-R1-core",
               "ROUND-1 primers — SPECIES-ID + RESISTANCE + DIAGNOSTIC/CNV (universal, all platforms)",
               "One-time order, reused every run. Pool-2 wells are shaded — run them as a separate small 3-plex.",
               R1_HDR, R1_W)
er, ncore = emit_r1(ws, sr, core, "PANEL_CORE")
ws.cell(er + 1, 3, f"Total oligos: {ncore}   ({len(core)} amplicons × 2)").font = Font(bold=True, color=INK)

ws, sr = sheet(wb, "Oligos-R1-diversity",
               "ROUND-1 primers — DIVERSITY module (COI / relatedness; 6 taxa)",
               "One-time order, reused every run. Co-pools with the core set into Pool 1 (single-pool panel).",
               R1_HDR, R1_W)
er, ndiv = emit_r1(ws, sr, div, "PANEL_DIV")
ws.cell(er + 1, 3, f"Total oligos: {ndiv}   ({len(div)} amplicons × 2)").font = Font(bold=True, color=INK)

# ======================================================================= R2 BARCODE SHEETS (combinatorial)
FWD_BCS = [bid for bid, _ in bcs[:12]]     # BC001..BC012  -> forward primers
REV_BCS = [bid for bid, _ in bcs[12:20]]   # BC013..BC020  -> reverse primers
R2_HDR = ["Plate", "Well", "Oligo_Name", "Sequence_5to3", "Length_nt", "Scale", "Purification", "Role", "Notes"]
R2_W   = [14, 7, 24, 78, 10, 10, 26, 10, 30]


def emit_r2(ws, start_row, platform):
    wells = wells_colmajor(); wi = 0; r = start_row
    rows = []
    if platform == "ONT":
        for bid in FWD_BCS:
            rows.append((f"R2-ONT-{bid}-F", bcmap[bid] + CS1, "forward (bc+CS1)",
                         "one-time; combine with any reverse"))
        for bid in REV_BCS:
            rows.append((f"R2-ONT-{bid}-R", bcmap[bid] + CS2, "reverse (bc+CS2)",
                         "one-time; combine with any forward"))
        purif = "Standard desalting (HPLC optional)"
    else:  # Illumina
        for bid in FWD_BCS:
            rows.append((f"R2-ILL-{bid}-i5", P5 + bcmap[bid] + CS1, "i5 / P5 (forward)",
                         "one-time; combine with any i7"))
        for bid in REV_BCS:
            rows.append((f"R2-ILL-{bid}-i7", P7 + bcmap[bid] + CS2, "i7 / P7 (reverse)",
                         "one-time; combine with any i5"))
        purif = "HPLC recommended (>60 nt)"
    for name, seq, role, note in rows:
        vals = [f"R2_{platform}_p1", wells[wi], name, seq, len(seq), "25 nmol", purif, role, note]
        for j, v in enumerate(vals, 1):
            cell = ws.cell(r, j, v)
            cell.font = f_mono if j in (2, 3, 4) else f_cell
            cell.alignment = wrap if j == 9 else top
            cell.border = border
        wi += 1; r += 1
    return r, len(rows)


ws, sr = sheet(wb, "Oligos-R2-barcode-ONT",
               "ROUND-2 barcoding primers — ONT (combinatorial dual index)",
               "12 forward (BC001–012) × 8 reverse (BC013–020) = 96 sample tags from 20 oligos. Reused every run.",
               R2_HDR, R2_W)
er, n = emit_r2(ws, sr, "ONT")
ws.cell(er + 1, 3, f"Total oligos: {n}  (12 fwd + 8 rev → 96 combinations)").font = Font(bold=True, color=INK)

ws, sr = sheet(wb, "Oligos-R2-index-Illumina",
               "ROUND-2 index primers — Illumina P5/P7 (combinatorial dual index)",
               "12 i5 (BC001–012) × 8 i7 (BC013–020) = 96 sample tags from 20 oligos. Reused every run.",
               R2_HDR, R2_W)
er, n = emit_r2(ws, sr, "ILL")
ws.cell(er + 1, 3, f"Total oligos: {n}  (12 i5 + 8 i7 → 96 combinations)").font = Font(bold=True, color=INK)

# ======================================================================= BARCODE LAYOUT / DEMUX MAP
ws, sr = sheet(wb, "Barcode-layout-96",
               "SAMPLE → BARCODE LAYOUT (96 samples) — also the demux --barcode_map input",
               "Each sample = one (fwd_bc, rev_bc) combination. Columns sample,fwd_bc,rev_bc are the exact demux map.",
               ["sample_plate_well", "sample", "fwd_bc", "rev_bc", "fwd_seq", "rev_seq"],
               [16, 16, 10, 10, 28, 28])
wells = wells_colmajor(); r = sr; s = 1
for fwd in FWD_BCS:            # 12
    for rev in REV_BCS:       # 8  -> 96
        w = wells[s - 1]
        vals = [w, f"S{s:03d}", fwd, rev, bcmap[fwd], bcmap[rev]]
        for j, v in enumerate(vals, 1):
            cell = ws.cell(r, j, v)
            cell.font = f_mono if j in (1, 3, 4, 5, 6) else f_cell
            cell.border = border
        r += 1; s += 1
ws.cell(r + 1, 2, "96 samples (12 fwd × 8 rev). Rename S001… to your sample IDs before ordering/demux.").font = f_note

# ======================================================================= REAGENT SHEETS
RG_HDR = ["Item", "Vendor", "Catalog # (indicative)", "Vendor URL", "Qty / notes (per 96-sample run)", "Verify"]
RG_W   = [46, 24, 26, 26, 46, 10]


def reagent_sheet(name, title, subtitle, rows):
    ws, r = sheet(wb, name, title, subtitle, RG_HDR, RG_W)
    for row in rows:
        note = row[0].startswith("NOTE")
        for j, v in enumerate(row, 1):
            cell = ws.cell(r, j, v)
            cell.font = f_note if note else (f_mono if j == 4 else f_cell)
            cell.alignment = wrap; cell.border = border
        r += 1
    return ws


reagent_sheet("Reagents-Common",
    "COMMON reagents & consumables (all platforms)",
    "Sample prep → two-step tailed-amplicon library. Quantities sized for 96 samples × 2 pools (~200 R1 rxns).",
    [
     ["Whatman 903 Protein Saver cards", "Cytiva", "10534612", "cytivalifesciences.com", "DBS collection (96 samples)", "confirm"],
     ["Harris Uni-Core 3.0 mm DBS punch", "Cytiva/GE", "WB100060", "cytivalifesciences.com", "DBS punches", "confirm"],
     ["Chelex 100 Resin, mol. biol. grade", "Bio-Rad", "1421253", "bio-rad.com", "5-min quick lysis", "confirm"],
     ["Saponin", "Sigma-Aldrich (Merck)", "47036", "sigmaaldrich.com", "quick-lysis (RBC lysis)", "confirm"],
     ["QIAamp DNA Mini Kit (optional ref. extraction)", "QIAGEN", "51304", "qiagen.com", "reference arm vs quick-lysis", "confirm"],
     ["QIAGEN Multiplex PCR Plus Kit", "QIAGEN", "206152", "qiagen.com", "Round-1 multiplex master mix; ~200 rxns → 2×(100 rxn) or 1×(1000 rxn, 206153)", "confirm"],
     ["Lyo-Ready Hot Start PCR mix (cold-chain-free alt.)", "Meridian Bioscience", "(contact vendor)", "meridianbioscience.com", "lyophilizable route; validate vs wet", "confirm"],
     ["NEB Q5 Hot Start HiFi 2X Master Mix", "New England Biolabs", "M0494", "neb.com", "Round-2 barcoding PCR; ~110 rxns → 1×(100 rxn, 5 mL)", "confirm"],
     ["AMPure XP beads (60 mL)", "Beckman Coulter", "A63881", "beckman.com", "SPRI clean-up / size-select (covers many runs)", "confirm"],
     ["Qubit dsDNA HS Assay Kit (500 assays)", "Thermo Fisher", "Q32854", "thermofisher.com", "quantification", "confirm"],
     ["Qubit 4 Fluorometer", "Thermo Fisher", "Q33238", "thermofisher.com", "one-time instrument", "confirm"],
     ["TapeStation / Bioanalyzer HS DNA (optional QC)", "Agilent", "5067-5584", "agilent.com", "amplicon-size QC of pools", "confirm"],
     ["Nuclease-free water", "Thermo Fisher", "AM9937", "thermofisher.com", "", "confirm"],
     ["96-well PCR plates + adhesive/foil seals", "various", "-", "-", "≥4 plates (R1 pool1, R1 pool2, R2, dilutions)", "confirm"],
     ["Filter tips / low-bind tubes", "various", "-", "-", "consumables", "confirm"],
    ])

reagent_sheet("Reagents-ONT",
    "ONT platform reagents (primary / clinic)",
    "R10.4.1 chemistry. One MinION flow cell comfortably runs a 96-plex amplicon pool.",
    [
     ["MinION Mk1B/Mk1C (or GridION)", "Oxford Nanopore", "(device)", "nanoporetech.com", "one-time instrument", "confirm"],
     ["MinION Flow Cell R10.4.1", "Oxford Nanopore", "FLO-MIN114", "nanoporetech.com", "1 per 96-sample run (buy 1–2)", "confirm"],
     ["Flongle Flow Cell + adapter (small runs)", "Oxford Nanopore", "FLO-FLG114 / ADP-FLG001", "nanoporetech.com", "cheaper for <24-sample runs", "confirm"],
     ["Ligation Sequencing Kit V14", "Oxford Nanopore", "SQK-LSK114", "nanoporetech.com", "adapter (motor) + reagents; ~ up to 6 preps/kit", "confirm"],
     ["NEBNext Companion Module for ONT", "New England Biolabs", "E7180S", "neb.com", "end-prep + adapter ligation", "confirm"],
     ["Flow Cell Wash Kit", "Oxford Nanopore", "EXP-WSH004", "nanoporetech.com", "reuse flow cells across runs", "confirm"],
     ["NOTE: barcodes are our own 24-nt oligos (Oligos-R2-barcode-ONT) — ONT PCR-barcoding kit NOT required.", "", "", "", "", "-"],
    ])

reagent_sheet("Reagents-Illumina",
    "Illumina / short-read reagents (reference lab)",
    "No library-prep kit — P5/P7 + index are added by the Round-2 index PCR.",
    [
     ["MiSeq / iSeq 100 / MiniSeq", "Illumina", "(device)", "illumina.com", "one-time instrument", "confirm"],
     ["MiSeq Reagent Kit v2 (500-cycle, 2×250)", "Illumina", "MS-102-2003", "illumina.com", "1 per run (v3 600-cyc MS-102-3003 for 2×300)", "confirm"],
     ["PhiX Control v3", "Illumina", "FC-110-3001", "illumina.com", "spike-in for low-diversity amplicons", "confirm"],
     ["KAPA Library Quantification Kit (Illumina)", "Roche", "KK4824", "sequencing.roche.com", "optional qPCR quant before loading", "confirm"],
     ["NOTE: NO library-prep kit needed — indices added by round-2 index PCR (Oligos-R2-index-Illumina).", "", "", "", "", "-"],
     ["NOTE: use CS1/CS2 as custom Read1/Read2/index primers (spike into cartridge), OR swap tails to TruSeq for stock primers.", "", "", "", "", "-"],
    ])

reagent_sheet("Reagents-MGI",
    "MGI / BGI platform reagents (reference lab)",
    "Round-1 amplicons are identical; MGI needs an adapter-conversion + circularization step.",
    [
     ["DNBSEQ-G99 / G400", "MGI Tech", "(device)", "en.mgi-tech.com", "one-time instrument", "confirm"],
     ["DNBSEQ-G99RS Sequencing Kit (FCL / App-A)", "MGI Tech", "(chemistry-dependent)", "en.mgi-tech.com", "read-length/kit per run", "confirm"],
     ["MGIEasy Circularization Kit", "MGI Tech", "1000005259", "en.mgi-tech.com", "convert PCR amplicons to DNBs", "confirm"],
     ["Make-DNB / DNB loading reagents", "MGI Tech", "(kit-dependent)", "en.mgi-tech.com", "DNB prep + loading", "confirm"],
     ["NOTE: convert via MGIEasy/App-A, OR design MGI-adapter round-2 primers (contact MGI for adapter sequences).", "", "", "", "", "-"],
    ])

# ======================================================================= CONTROLS & STRAINS
ws, r = sheet(wb, "Controls-and-strains",
    "CONTROLS, STRAINS & CHARACTERIZED SAMPLES (6 taxa)",
    "Cultured strains where they exist; WHO/NIBSC DNA standards + synthetic gBlocks for the non-culturable species and defined non-Pf resistance alleles.",
    ["Item", "Source", "Catalog # (indicative)", "URL", "Role in run", "Verify"],
    [46, 26, 26, 24, 50, 10])
controls = [
 ["P. falciparum 3D7", "BEI Resources (MR4)", "MRA-102", "beiresources.org", "CNV 1-copy calibrator; hrp2+/hrp3+; dhps A437G; WT pfcrt", "confirm"],
 ["P. falciparum Dd2", "BEI Resources (MR4)", "MRA-150", "beiresources.org", "pfmdr1 ×4 (CNV+); pfhrp2 DELETED; pfcrt CVIET", "confirm"],
 ["P. falciparum HB3", "BEI Resources (MR4)", "MRA-155", "beiresources.org", "pfhrp3 DELETED control", "confirm"],
 ["P. falciparum 7G8", "BEI Resources (MR4)", "MRA-152", "beiresources.org", "pfcrt SVMNT haplotype", "confirm"],
 ["P. falciparum D10", "BEI Resources (MR4)", "MRA-201", "beiresources.org", "2nd pfhrp2-deleted control", "confirm"],
 ["P. vivax genomic DNA (Sal-1 / Pacific)", "BEI Resources / NIBSC", "(confirm)", "beiresources.org", "Pv species + pvmdr1/pvdhfr/pvdhps positive control", "confirm"],
 ["P. knowlesi A1-H.1 gDNA (culture-adapted)", "collaborators / EBI / NIBSC", "(confirm)", "-", "Pk species-ID + Pk cluster typing positive control", "confirm"],
 ["WHO Int'l Standard — P. falciparum DNA", "NIBSC", "04/176", "nibsc.org", "quantitative Pf reference", "confirm"],
 ["WHO/NIBSC malaria DNA reference panel (Pv/Pm/Po/Pk)", "NIBSC", "(confirm codes)", "nibsc.org", "species-ID positive controls for non-culturable taxa", "confirm"],
 ["Synthetic control — cox3 species-ID, 6 taxa (gBlocks)", "IDT", "(custom)", "idtdna.com", "defined 6-taxa species-ID positive control", "-"],
 ["Synthetic control — non-Pf resistance alleles (gBlocks)", "IDT", "(custom)", "idtdna.com", "pmdhfr 57/58/114, podhfr 15/57/58/113 + discovery dhps/crt/mdr1/k13 windows (incl. PmUG01 polarity S vs R/N)", "-"],
 ["Synthetic control — Pk cluster diagnostic SNPs (gBlocks)", "IDT", "(custom)", "idtdna.com", "Mf / Mn / cluster-3 diagnostic alleles", "-"],
 ["Characterized field DBS (per species, where available)", "local biobank / collaborators", "-", "-", "real-sample arm", "-"],
 ["No-template control (NTC)", "-", "-", "-", "contamination monitor (assign a barcode in Barcode-layout-96)", "-"],
]
for row in controls:
    for j, v in enumerate(row, 1):
        cell = ws.cell(r, j, v)
        cell.font = f_mono if j == 4 else f_cell
        cell.alignment = wrap; cell.border = border
    r += 1

# ---- order sheets ----
order = ["READ ME", "Oligos-R1-core", "Oligos-R1-diversity", "Oligos-R2-barcode-ONT",
         "Oligos-R2-index-Illumina", "Barcode-layout-96", "Reagents-Common", "Reagents-ONT",
         "Reagents-Illumina", "Reagents-MGI", "Controls-and-strains"]
wb._sheets.sort(key=lambda s: order.index(s.title))

out = f"{HERE}/PANPLASMO-ORDER-WORKBOOK.xlsx"
wb.save(out)
print("wrote", out)
print("sheets:", [s.title for s in wb.worksheets])
print(f"R1 core oligos={ncore}, R1 diversity oligos={ndiv}, total R1={ncore+ndiv}")
