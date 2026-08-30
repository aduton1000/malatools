#!/usr/bin/env python3
"""
add_cost_estimate.py — append a 'Cost-estimate' sheet sizing each platform's run to its optimum multiplexing
for the 289-amplicon panel, and cost one-time vs per-run buckets. Indicative USD (economy oligo/kit tier),
every figure editable + flagged VERIFY. Re-runnable (drops an existing Cost-estimate sheet first).
"""
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
import os

HERE = os.path.dirname(os.path.abspath(__file__))
WB = f"{HERE}/PANPLASMO-ORDER-WORKBOOK.xlsx"

# ---------- panel + depth model ----------
N_AMP      = 289
AMP_BP     = 300               # mean informative length
TARGET_X   = 400               # comfortable mean depth per amplicon (floor ~100x on smallest)
BASES_PER_SAMPLE = N_AMP * AMP_BP * TARGET_X          # ~34.7 Mb
READS_PER_SAMPLE = N_AMP * TARGET_X                   # ~115.6k merged reads

# ---------- platform yields (usable, conservative, indicative) ----------
# (label, usable_bases, usable_reads_or_None, chosen_optimum_samples, barcode_oligos_needed, notes)
PLATFORMS = [
    ("ONT MinION R10.4.1", 15e9, 30e6, 384, 40,
     "≈15 Gb usable → depth-limited ~428; set to 384 (4 plates) ≈450×. Push to 576 (~300×) for lowest $/sample."),
    ("ONT Flongle",         1.0e9, 2.5e6, 24, 0,
     "small/urgent batches; reuse a subset of the 96-barcode set (no extra oligos)."),
    ("Illumina MiSeq v2 (2×250)", 7.5e9, 15e6, 96, 20,
     "read-limited (~15M PE) → ~120; set to 96 (1 plate) ≈540×. MiSeq v3 600-cyc → ~150."),
    ("MGI DNBSEQ-G99 (PE)", 10e9, None, 192, 28,
     "FCL kit ≈10 Gb → ~285; set to 192 (16×12 barcodes) ≈430×. Scales with FCS/FCL choice → ~288."),
]

# ---------- unit prices (indicative USD, economy tier — EDIT/VERIFY) ----------
PX = {
    "oligo_desalt_per_base": 0.13,     # 25 nmol standard desalting, plate
    "oligo_hplc_per_base":   0.40,     # >60 nt HPLC (Illumina index)
    "prep_per_sample":       3.00,     # lysis + R1(2 pools) + R2 + SPRI + Qubit + plastics (bulk)
    # per-run sequencing (fixed, at optimum):
    "seq_ONT MinION R10.4.1": 1100,    # FLO-MIN114 ~900 + LSK114 share ~120 + companion/wash ~80
    "seq_ONT Flongle":        200,     # FLO-FLG114 ~120 + prep share
    "seq_Illumina MiSeq v2 (2×250)": 1440,  # MS-102-2003 ~1400 + PhiX ~40
    "seq_MGI DNBSEQ-G99 (PE)": 1250,   # FCL kit ~800 + circularization ~300 + DNB ~150
}
# one-time oligo bases (from the workbook)
wb = openpyxl.load_workbook(WB)
def bases(sheet):
    ws = wb[sheet]; b = 0
    for row in ws.iter_rows(min_row=5, values_only=True):
        if not row[2] or str(row[2]).startswith("Total"): continue
        b += int(row[4])
    return b
R1_BASES = bases("Oligos-R1-core") + bases("Oligos-R1-diversity")   # 25,895
R1_COST  = R1_BASES * PX["oligo_desalt_per_base"]
# barcode oligo cost per platform (ONT/MGI desalted 46 nt; Illumina HPLC ~73 nt)
def bc_cost(name, n_oligos):
    if n_oligos == 0:
        return 0.0
    if "Illumina" in name:
        return n_oligos * 73 * PX["oligo_hplc_per_base"]
    return n_oligos * 46 * PX["oligo_desalt_per_base"]

CONTROLS = 2505   # Pf strains(5×85) + Pv gDNA(250) + Pk gDNA(300) + NIBSC Pf 04/176(300) + NIBSC Pv/Pm/Po/Pk(600) + 3 gBlocks(450) + QIAamp(180)

# ---------- styling ----------
INK="231D19"; PAPER="FBF8F2"; BLOOD="9E2B25"; TEAL="0E6E66"; GOLD="8F6B1C"; BAND="EFEAE0"
f_title=Font(size=13,bold=True,color=BLOOD); f_sub=Font(size=10,italic=True,color="6C625A")
f_hdr=Font(size=10,bold=True,color=PAPER); f_cell=Font(size=10,color=INK)
f_b=Font(size=10,bold=True,color=INK); f_money=Font(size=10,color=INK); f_note=Font(size=9,italic=True,color=GOLD)
fill_hdr=PatternFill("solid",fgColor=INK); fill_band=PatternFill("solid",fgColor=BAND)
fill_tot=PatternFill("solid",fgColor="F6E7C9"); fill_teal=PatternFill("solid",fgColor="D8ECE9")
thin=Side(style="thin",color="DED4C5"); bd=Border(bottom=thin); wrap=Alignment(wrap_text=True,vertical="top")

if "Cost-estimate" in wb.sheetnames:
    del wb["Cost-estimate"]
ws = wb.create_sheet("Cost-estimate")
for j,w in enumerate([34,20,20,20,20,44],1):
    ws.column_dimensions[get_column_letter(j)].width=w

r=1
ws.cell(r,1,"COST ESTIMATE — run sized to each platform's optimum multiplexing").font=f_title; r+=1
ws.cell(r,1,f"Panel {N_AMP} amplicons · target ≥{TARGET_X}× mean/amplicon (~{BASES_PER_SAMPLE/1e6:.0f} Mb, "
            f"~{READS_PER_SAMPLE/1000:.0f}k reads per sample). Indicative USD, economy tier — EDIT prices below, VERIFY all.").font=f_sub
r+=2

# ---- optimum sizing table ----
ws.cell(r,1,"1 · OPTIMUM RUN SIZE PER PLATFORM").font=f_b; r+=1
hdr=["Platform","Usable yield","Samples / run (optimum)","Depth @ optimum","Barcode oligos","Basis / notes"]
for j,h in enumerate(hdr,1):
    c=ws.cell(r,j,h); c.font=f_hdr; c.fill=fill_hdr; c.alignment=wrap
r+=1
opt={}
for name,ybases,yreads,nsamp,nbc,note in PLATFORMS:
    depth=ybases/(nsamp*N_AMP*AMP_BP)
    opt[name]=(nsamp,nbc)
    vals=[name, f"{ybases/1e9:.1f} Gb", f"{nsamp}", f"~{depth:.0f}× mean", f"{nbc or '—'}", note]
    for j,v in enumerate(vals,1):
        c=ws.cell(r,j,v); c.font=f_cell if j!=3 else f_b; c.alignment=wrap; c.border=bd
    r+=1
r+=1

# ---- unit prices (editable) ----
ws.cell(r,1,"2 · UNIT PRICES (indicative USD — edit & verify)").font=f_b; r+=1
price_rows=[
    ("Oligo — 25 nmol desalted, per base", PX["oligo_desalt_per_base"]),
    ("Oligo — >60 nt HPLC (Illumina index), per base", PX["oligo_hplc_per_base"]),
    ("Library prep per sample (lysis+2×PCR+SPRI+Qubit+plastics)", PX["prep_per_sample"]),
    ("Controls & strains (one-time, semi-reusable)", CONTROLS),
]
for lbl,val in price_rows:
    ws.cell(r,1,lbl).font=f_cell; c=ws.cell(r,2,val); c.font=f_money; c.number_format='$#,##0.00'; r+=1
r+=1

# ---- one-time bucket ----
ws.cell(r,1,"3 · ONE-TIME LAB SETUP (buy once, reused every run)").font=f_b; ws.cell(r,1).fill=fill_band; r+=1
for j,h in enumerate(["Item","","ONT MinION","Illumina MiSeq v2","MGI DNBSEQ-G99",""],1):
    c=ws.cell(r,j,h); c.font=f_hdr; c.fill=fill_hdr
r+=1
onetime={}
def money(cell,v):
    cell.value=v; cell.number_format='$#,##0'; cell.font=f_money
# R1 shared
ws.cell(r,1,f"Round-1 primer panel (578 oligos, {R1_BASES:,} nt)").font=f_cell
for j in (3,4,5): money(ws.cell(r,j),R1_COST)
r+=1
# R2 barcodes per platform
bc_ont=bc_cost("ONT",opt["ONT MinION R10.4.1"][1])
bc_ill=bc_cost("Illumina",opt["Illumina MiSeq v2 (2×250)"][1])
bc_mgi=bc_cost("MGI",opt["MGI DNBSEQ-G99 (PE)"][1])
ws.cell(r,1,"Round-2 barcode/index oligos (sized to optimum)").font=f_cell
money(ws.cell(r,3),bc_ont); money(ws.cell(r,4),bc_ill); money(ws.cell(r,5),bc_mgi); r+=1
ws.cell(r,1,"Controls, strains & gBlocks").font=f_cell
for j in (3,4,5): money(ws.cell(r,j),CONTROLS)
r+=1
# subtotal
ws.cell(r,1,"One-time subtotal").font=f_b
ot_ont=R1_COST+bc_ont+CONTROLS; ot_ill=R1_COST+bc_ill+CONTROLS; ot_mgi=R1_COST+bc_mgi+CONTROLS
onetime={"ONT":ot_ont,"ILL":ot_ill,"MGI":ot_mgi}
for j,v in zip((3,4,5),(ot_ont,ot_ill,ot_mgi)):
    money(ws.cell(r,j),v); ws.cell(r,j).font=f_b; ws.cell(r,j).fill=fill_tot
r+=2

# ---- per-run bucket ----
ws.cell(r,1,"4 · PER-RUN (recurring, at platform optimum)").font=f_b; ws.cell(r,1).fill=fill_band; r+=1
for j,h in enumerate(["Item","","ONT MinION","Illumina MiSeq v2","MGI DNBSEQ-G99",""],1):
    c=ws.cell(r,j,h); c.font=f_hdr; c.fill=fill_hdr
r+=1
n_ont=opt["ONT MinION R10.4.1"][0]; n_ill=opt["Illumina MiSeq v2 (2×250)"][0]; n_mgi=opt["MGI DNBSEQ-G99 (PE)"][0]
ws.cell(r,1,"Samples this run").font=f_cell
for j,v in zip((3,4,5),(n_ont,n_ill,n_mgi)):
    ws.cell(r,j,v).font=f_b
r+=1
prep_ont=n_ont*PX["prep_per_sample"]; prep_ill=n_ill*PX["prep_per_sample"]; prep_mgi=n_mgi*PX["prep_per_sample"]
ws.cell(r,1,"Library prep (all samples)").font=f_cell
money(ws.cell(r,3),prep_ont); money(ws.cell(r,4),prep_ill); money(ws.cell(r,5),prep_mgi); r+=1
seq_ont=PX["seq_ONT MinION R10.4.1"]; seq_ill=PX["seq_Illumina MiSeq v2 (2×250)"]; seq_mgi=PX["seq_MGI DNBSEQ-G99 (PE)"]
ws.cell(r,1,"Sequencing (flow cell / kit)").font=f_cell
money(ws.cell(r,3),seq_ont); money(ws.cell(r,4),seq_ill); money(ws.cell(r,5),seq_mgi); r+=1
run_ont=prep_ont+seq_ont; run_ill=prep_ill+seq_ill; run_mgi=prep_mgi+seq_mgi
ws.cell(r,1,"Per-run subtotal").font=f_b
for j,v in zip((3,4,5),(run_ont,run_ill,run_mgi)):
    money(ws.cell(r,j),v); ws.cell(r,j).font=f_b; ws.cell(r,j).fill=fill_tot
r+=1
ws.cell(r,1,"→ Steady-state cost / sample").font=f_b
for j,v,n in zip((3,4,5),(run_ont,run_ill,run_mgi),(n_ont,n_ill,n_mgi)):
    c=ws.cell(r,j,v/n); c.number_format='$#,##0.00'; c.font=f_b; c.fill=fill_teal
r+=2

# ---- headline totals ----
ws.cell(r,1,"5 · WHAT YOU PAY").font=f_b; ws.cell(r,1).fill=fill_band; r+=1
ws.cell(r,1,"FIRST run all-in (one-time + one optimum run; instrument owned)").font=f_cell
for j,ot,run in zip((3,4,5),(ot_ont,ot_ill,ot_mgi),(run_ont,run_ill,run_mgi)):
    money(ws.cell(r,j),ot+run); ws.cell(r,j).font=f_b
r+=1
ws.cell(r,1,"EVERY later run (oligos+controls already owned)").font=f_cell
for j,run,n in zip((3,4,5),(run_ont,run_ill,run_mgi),(n_ont,n_ill,n_mgi)):
    money(ws.cell(r,j),run);
r+=1
ws.cell(r,1,"Effective $/sample after setup").font=f_cell
for j,run,n in zip((3,4,5),(run_ont,run_ill,run_mgi),(n_ont,n_ill,n_mgi)):
    c=ws.cell(r,j,run/n); c.number_format='$#,##0.00'
r+=2

notes=[
 "Instruments (capital) are EXCLUDED — MinION ~$1,000; MiSeq/DNBSEQ-G99 are lab instruments. Add if not owned.",
 "One-time oligos+controls amortize to near-zero per sample across many runs (e.g. ONT $%.0f over 10×384 = $%.2f/sample)."%(ot_ont, ot_ont/(10*n_ont)),
 "MiSeq v2 is the pricey short-read option at its optimum (fixed cartridge over 96) — move to NextSeq/NovaSeq at their optimum (100s–1000s of samples) to reach ~$2–4/sample.",
 "ONT and MGI at optimum meet the $1–10/sample target; ONT is cheapest at small-to-mid scale, MGI/short-read cheapest at very large scale.",
 "All prices indicative (2026, economy tier); confirm against current quotes. IDT/premium oligos run ~1.8–2× the economy per-base used here.",
]
for n in notes:
    c=ws.cell(r,1,"• "+n); c.font=f_note; c.alignment=wrap; ws.merge_cells(start_row=r,start_column=1,end_row=r,end_column=6); r+=1

# move Cost-estimate to position 2 (after READ ME)
order=[s.title for s in wb.worksheets if s.title not in ("Cost-estimate",)]
order.insert(1,"Cost-estimate")
wb._sheets.sort(key=lambda s: order.index(s.title))
wb.save(WB)

# ---- console summary ----
print("OPTIMUM SIZING & COST (indicative USD)")
for tag,n,ot,run in [("ONT MinION",n_ont,ot_ont,run_ont),("Illumina MiSeq v2",n_ill,ot_ill,run_ill),("MGI DNBSEQ-G99",n_mgi,ot_mgi,run_mgi)]:
    print(f"  {tag:20} {n:>4} samples/run | one-time ${ot:>6,.0f} | per-run ${run:>6,.0f} | ${run/n:5.2f}/sample | first-run all-in ${ot+run:>6,.0f}")
print(f"  R1 panel one-time: ${R1_COST:,.0f}  | controls one-time: ${CONTROLS:,.0f}")
