#!/usr/bin/env python3
"""add_cost_summary.py — append a clean 'Cost-summary' sheet: 3 platforms side by side at their optimum
multiplexing, best value per row highlighted, plus a decision guide. Re-runnable. Figures mirror
add_cost_estimate.py exactly."""
import openpyxl, os
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

HERE = os.path.dirname(os.path.abspath(__file__)); WB = f"{HERE}/PANPLASMO-ORDER-WORKBOOK.xlsx"
N_AMP=289; AMP_BP=300; TARGET_X=400
PX={"desalt":0.13,"hplc":0.40,"prep":3.00}
CONTROLS=2505
wb=openpyxl.load_workbook(WB)
def bases(sheet):
    ws=wb[sheet]; b=0
    for row in ws.iter_rows(min_row=5,values_only=True):
        if not row[2] or str(row[2]).startswith("Total"): continue
        b+=int(row[4])
    return b
R1=(bases("Oligos-R1-core")+bases("Oligos-R1-diversity"))*PX["desalt"]

# platform: name, samples, barcode_oligos, bc_hplc?, seq_fixed
P=[("ONT MinION",384,40,False,1100),
   ("Illumina MiSeq v2",96,20,True,1440),
   ("MGI DNBSEQ-G99",192,28,False,1250)]
rows=[]
for name,n,nbc,hplc,seq in P:
    bc=nbc*(73*PX["hplc"] if hplc else 46*PX["desalt"])
    onetime=R1+bc+CONTROLS
    prep=n*PX["prep"]; run=prep+seq
    rows.append(dict(name=name,n=n,onetime=onetime,run=run,persample=run/n,first=onetime+run,
                     depth=15e9/(n*N_AMP*AMP_BP) if name.startswith("ONT") else (7.5e9 if "MiSeq" in name else 10e9)/(n*N_AMP*AMP_BP)))

# styling
INK="231D19";PAPER="FBF8F2";BLOOD="9E2B25";TEAL="0E6E66";GOLD="8F6B1C";BAND="EFEAE0";GREEN="1E7A3D"
f_title=Font(size=13,bold=True,color=BLOOD);f_sub=Font(size=10,italic=True,color="6C625A")
f_hdr=Font(size=11,bold=True,color=PAPER);f_metric=Font(size=10,bold=True,color=INK)
f_cell=Font(size=11,color=INK);f_best=Font(size=11,bold=True,color=GREEN);f_note=Font(size=9,italic=True,color=GOLD)
fill_hdr=PatternFill("solid",fgColor=INK);fill_best=PatternFill("solid",fgColor="DDF0E1")
fill_band=PatternFill("solid",fgColor=BAND);fill_rec=PatternFill("solid",fgColor="F6E7C9")
thin=Side(style="thin",color="DED4C5");bd=Border(bottom=thin,left=thin,right=thin,top=thin)
ctr=Alignment(horizontal="center",vertical="center",wrap_text=True);left=Alignment(vertical="center",wrap_text=True)

if "Cost-summary" in wb.sheetnames: del wb["Cost-summary"]
ws=wb.create_sheet("Cost-summary")
ws.column_dimensions["A"].width=30
for col in "BCD": ws.column_dimensions[col].width=22
ws.column_dimensions["E"].width=2

r=1
ws.cell(r,1,"COST SUMMARY — three platforms compared").font=f_title; r+=1
ws.cell(r,1,f"Each platform at its OPTIMUM multiplexing for the {N_AMP}-amplicon panel (≥{TARGET_X}× mean depth). "
            "Indicative USD; green = best in row. Full breakdown in 'Cost-estimate'.").font=f_sub
ws.merge_cells(start_row=r,start_column=1,end_row=r,end_column=4); r+=2

# header row
hdr=["",rows[0]["name"],rows[1]["name"],rows[2]["name"]]
for j,h in enumerate(hdr,1):
    c=ws.cell(r,j,h); c.font=f_hdr; c.fill=fill_hdr; c.alignment=ctr; c.border=bd
ws.cell(r,1,"Metric").font=f_hdr; ws.cell(r,1).fill=fill_hdr; ws.cell(r,1).alignment=ctr
r+=1

def band(metric,key,fmt,best="min",money=True,intval=False):
    global r
    vals=[x[key] for x in rows]
    pick=min(vals) if best=="min" else max(vals)
    ws.cell(r,1,metric).font=f_metric; ws.cell(r,1).alignment=left; ws.cell(r,1).border=bd
    for j,x in enumerate(vals,2):
        c=ws.cell(r,j)
        if intval: c.value=int(x)
        else: c.value=x
        c.number_format=fmt; c.alignment=ctr; c.border=bd
        if abs(x-pick)<1e-6:
            c.font=f_best; c.fill=fill_best
        else:
            c.font=f_cell
    r+=1

band("Samples per run (optimum)","n","#,##0",best="max",intval=True)
band("Mean depth / amplicon","depth","#,##0\"×\"",best="max")
band("One-time setup (oligos+controls)","onetime","$#,##0",best="min")
band("Per-run (prep + sequencing)","run","$#,##0",best="min")
band("Cost per sample (steady-state)","persample","$#,##0.00",best="min")
band("First run all-in (setup + 1 run)","first","$#,##0",best="min")
r+=1

# bottom line per platform
ws.cell(r,1,"BOTTOM LINE").font=f_metric; ws.cell(r,1).fill=fill_band
ws.merge_cells(start_row=r,start_column=1,end_row=r,end_column=4); r+=1
bl=["Cheapest per sample & lowest entry cost — best for routine surveillance and small-to-mid batches. Same-day turnaround.",
    "Highest per-sample at its MiSeq optimum (fixed cartridge / 96). Use only if a MiSeq is the on-hand instrument; scale to NextSeq/NovaSeq for cheap short-read.",
    "Strong mid-scale short-read option; meets the $1–10 target at 192/run. Needs the extra circularization step."]
for x,txt in zip(rows,bl):
    ws.cell(r,1,x["name"]).font=Font(size=10,bold=True,color=TEAL); ws.cell(r,1).alignment=left
    c=ws.cell(r,2,txt); c.font=f_cell; c.alignment=left; ws.merge_cells(start_row=r,start_column=2,end_row=r,end_column=4)
    ws.row_dimensions[r].height=42
    r+=1
r+=1

# recommendation
ws.cell(r,1,"RECOMMENDATION").font=f_metric; ws.cell(r,1).fill=fill_rec
ws.merge_cells(start_row=r,start_column=1,end_row=r,end_column=4); r+=1
rec=[f"Primary / clinic: ONT MinION @384 → ${rows[0]['persample']:.2f}/sample (≈$4 at 576-plex). Fastest, cheapest to start.",
     f"Mid-scale reference lab: MGI DNBSEQ-G99 @192 → ${rows[2]['persample']:.2f}/sample.",
     f"Very large short-read batches: move off MiSeq to NextSeq/NovaSeq at its optimum → ~$2–4/sample.",
     "One-time oligos+controls (~$6k) amortize to near-zero per sample across runs; instruments excluded (add if not owned).",
     "All figures indicative (2026, economy tier) — verify against live quotes; premium oligos ~1.8–2× higher."]
for t in rec:
    c=ws.cell(r,1,"• "+t); c.font=f_note; c.alignment=left; ws.merge_cells(start_row=r,start_column=1,end_row=r,end_column=4); r+=1

# place Cost-summary right after READ ME
order=[s.title for s in wb.worksheets if s.title!="Cost-summary"]
order.insert(1,"Cost-summary")
wb._sheets.sort(key=lambda s: order.index(s.title))
wb.save(WB)
print("added Cost-summary. sheets:",[s.title for s in wb.worksheets])
for x in rows: print(f"  {x['name']:20} {x['n']:>4} | one-time ${x['onetime']:>6,.0f} | run ${x['run']:>6,.0f} | ${x['persample']:5.2f}/samp | first ${x['first']:>6,.0f}")
