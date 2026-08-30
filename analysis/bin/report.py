#!/usr/bin/env python3
"""
report.py — per-sample clinical HTML report from the analysis TSVs.
Usage: report.py <sample_outdir> <sample>
Writes <sample_outdir>/<sample>.report.html
"""
import sys, os, csv, html

PALETTE = """
:root{--paper:#F5F1EB;--surface:#FCFAF6;--ink:#231D19;--muted:#6C625A;--blood:#9E2B25;--teal:#0E6E66;
 --gold:#8F6B1C;--blue:#3B6EA5;--line:#DED4C5;}
*{box-sizing:border-box;margin:0;padding:0}
body{background:var(--paper);color:var(--ink);font-family:"IBM Plex Sans",system-ui,sans-serif;line-height:1.5;padding:24px}
.wrap{max-width:900px;margin:0 auto}
.ey{font:600 11px monospace;letter-spacing:.16em;color:var(--blood);text-transform:uppercase}
h1{font:700 26px Georgia,serif;margin:2px 0 2px}
.sub{color:var(--muted);font-size:13px;margin-bottom:16px}
h2{font:700 12px monospace;letter-spacing:.06em;text-transform:uppercase;color:var(--blood);
 border-bottom:1px solid var(--line);padding-bottom:5px;margin:20px 0 8px}
table{width:100%;border-collapse:collapse;font-size:12.5px;margin:6px 0}
th{background:var(--ink);color:#fff;text-align:left;padding:5px 8px;font:600 10px monospace;text-transform:uppercase}
td{padding:5px 8px;border-bottom:1px solid var(--line);vertical-align:top}
tr:nth-child(even) td{background:var(--surface)}
.pill{display:inline-block;font:600 10px monospace;padding:2px 8px;border-radius:10px}
.R{background:#f3d9d6;color:var(--blood)} .S{background:#d7ece9;color:var(--teal)}
.warn{background:#f6e8cf;color:var(--gold)} .mut{color:var(--blood);font-weight:700}
.sptag{display:inline-block;background:var(--blue);color:#fff;font:600 11px monospace;padding:2px 9px;border-radius:10px;margin-right:6px}
.mixed{background:var(--gold)}
.foot{margin-top:24px;font:500 10px monospace;color:var(--muted);border-top:1px solid var(--line);padding-top:10px}
.none{color:var(--muted);font-style:italic;font-size:12px}
"""


def rd(p):
    return list(csv.reader(open(p), delimiter="\t")) if os.path.exists(p) else []


def main():
    od, s = sys.argv[1], sys.argv[2]
    sp = rd(f"{od}/{s}.species.tsv")[1:]
    res = rd(f"{od}/{s}.resistance.tsv")[1:]
    drugs = rd(f"{od}/{s}.drugs.tsv")[1:]
    cnv = rd(f"{od}/{s}.cnv.tsv")[1:]
    div = rd(f"{od}/{s}.diversity.tsv")[1:]
    clus = rd(f"{od}/{s}.cluster.tsv")

    present = [(r[1], float(r[3])) for r in sp if r[4] == "present"]
    present.sort(key=lambda x: -x[1])
    mixed = len(present) > 1
    E = html.escape
    out = [f"<!doctype html><html><head><meta charset=utf-8><title>{E(s)} report</title>",
           '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600;700&display=swap">',
           f"<style>{PALETTE}</style></head><body><div class=wrap>",
           '<div class=ey>MalaCall · clinical report</div>',
           f"<h1>{E(s)}</h1>"]

    # species
    if present:
        tag = " ".join(f'<span class="sptag{" mixed" if mixed else ""}">{E(x[0])} {x[1]*100:.0f}%</span>' for x in present)
        out.append(f'<div class=sub>{"MIXED INFECTION · " if mixed else ""}{tag}</div>')
    else:
        out.append('<div class=sub><b>No malaria parasite detected.</b></div>')

    # resistance per species
    out.append("<h2>Drug resistance</h2>")
    if not drugs:
        out.append('<div class=none>No resistance markers assessed.</div>')
    for spname, _ in present:
        dr = [r for r in drugs if r[1] == spname and r[3] != "no resistance marker detected"]
        out.append(f"<h3 style='font:600 13px sans-serif;margin:10px 0 3px'>{E(spname)}</h3>")
        if dr:
            out.append("<table><tr><th>Drug</th><th>Status</th><th>Markers</th></tr>")
            for r in dr:
                out.append(f"<tr><td>{E(r[2])}</td><td><span class='pill R'>{E(r[3])}</span></td><td class=mut>{E(r[4])}</td></tr>")
            out.append("</table>")
        else:
            out.append('<div class=none>No resistance-associated mutations detected.</div>')
    # detailed resistant codons
    resR = [r for r in res if r[9] == "RESISTANT"]
    if resR:
        out.append("<details><summary style='cursor:pointer;font-size:12px;color:var(--muted)'>Resistant codon detail</summary><table>"
                   "<tr><th>Species</th><th>Gene</th><th>Codon</th><th>Sensitive</th><th>Observed</th><th>VAF</th><th>Depth</th></tr>")
        for r in resR:
            out.append(f"<tr><td>{E(r[1])}</td><td>{E(r[2])}</td><td>{E(r[3])}</td><td>{E(r[5])}</td>"
                       f"<td class=mut>{E(r[6])}</td><td>{E(r[7])}</td><td>{E(r[8])}</td></tr>")
        out.append("</table></details>")

    # diagnostic / CNV
    out.append("<h2>Diagnostic resistance &amp; copy number</h2>")
    flags = [r for r in cnv if "DELETED" in r[7] or "AMPLIFIED" in r[7]]
    if flags:
        out.append("<table><tr><th>Species</th><th>Gene</th><th>Type</th><th>Call</th><th>Ratio</th><th>Significance</th></tr>")
        for r in flags:
            cls = "warn"
            out.append(f"<tr><td>{E(r[1])}</td><td>{E(r[2])}</td><td>{E(r[3])}</td>"
                       f"<td><span class='pill {cls}'>{E(r[7])}</span></td><td>{E(r[6])}</td><td>{E(r[8])}</td></tr>")
        out.append("</table>")
    else:
        out.append('<div class=none>No gene deletions or amplifications detected.</div>')

    # diversity
    out.append("<h2>Parasite diversity</h2>")
    if div:
        out.append("<table><tr><th>Species</th><th>Loci typed</th><th>Polymorphic</th><th>Mean He</th><th>COI (complexity of infection)</th></tr>")
        for r in div:
            out.append(f"<tr><td>{E(r[1])}</td><td>{E(r[2])}</td><td>{E(r[3])}</td><td>{E(r[5])}</td><td><b>{E(r[6])}</b></td></tr>")
        out.append("</table>")
    else:
        out.append('<div class=none>Insufficient depth for diversity typing.</div>')

    # Pk cluster
    if clus:
        call = next((r for r in clus if len(r) > 1 and r[1] == "CALL"), None)
        if call:
            out.append("<h2>P. knowlesi subpopulation</h2>")
            out.append(f"<div style='font-size:14px'><span class='sptag'>{E(call[5])}</span></div>")

    out.append('<div class=foot>MalaCall analysis pipeline · in-silico validated · '
               'resistance called vs the SENSITIVE allele · NOT FOR CLINICAL USE (research pipeline; wet-lab validation pending)</div>')
    out.append("</div></body></html>")
    with open(f"{od}/{s}.report.html", "w") as fh:
        fh.write("\n".join(out))
    sys.stderr.write(f"[report] {s}: {od}/{s}.report.html\n")


if __name__ == "__main__":
    main()
