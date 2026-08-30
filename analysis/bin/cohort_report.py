#!/usr/bin/env python3
"""
cohort_report.py — aggregate all per-sample outputs into a surveillance summary (species distribution,
per-drug resistance frequencies, RDT-escape deletion frequencies, gene amplification, COI distribution).
Usage: cohort_report.py <out_dir> <cohort_name>
Writes <out_dir>/<cohort_name>.cohort.tsv and .cohort.html
"""
import sys, os, csv, glob, html
from collections import defaultdict, Counter


def rd(p):
    return list(csv.reader(open(p), delimiter="\t"))[1:] if os.path.exists(p) else []


def main():
    out, name = sys.argv[1], sys.argv[2]
    samples = sorted(d for d in os.listdir(out) if os.path.isdir(os.path.join(out, d)))
    n = len(samples)
    sp_count = Counter(); drug = defaultdict(lambda: defaultdict(int)); drug_tot = defaultdict(lambda: defaultdict(int))
    deletion = defaultdict(lambda: [0, 0]); cnv = defaultdict(lambda: [0, 0]); coi = defaultdict(list)
    for s in samples:
        od = os.path.join(out, s)
        present = [r[1] for r in rd(f"{od}/{s}.species.tsv") if len(r) > 4 and r[4] == "present"]
        for sp in present:
            sp_count[sp] += 1
        for r in rd(f"{od}/{s}.drugs.tsv"):
            sp, dg, status = r[1], r[2], r[3]
            if sp in present:
                drug_tot[sp][dg] += 1
                if "RESISTANCE" in status:
                    drug[sp][dg] += 1
        for r in rd(f"{od}/{s}.cnv.tsv"):
            sp, gene, typ, call = r[1], r[2], r[3], r[7]
            if sp not in present:
                continue
            if typ == "del":
                deletion[gene][1] += 1
                if "DELETED" in call:
                    deletion[gene][0] += 1
            elif typ == "cnv":
                cnv[gene][1] += 1
                if "AMPLIFIED" in call:
                    cnv[gene][0] += 1
        for r in rd(f"{od}/{s}.diversity.tsv"):
            coi[r[1]].append(int(r[6]))

    # TSV
    with open(f"{out}/{name}.cohort.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "key", "n_positive", "n_assessed", "frequency"])
        for sp, c in sp_count.most_common():
            w.writerow(["species", sp, c, n, round(c / n, 3)])
        for sp in drug:
            for dg, c in sorted(drug[sp].items()):
                t = drug_tot[sp][dg]
                w.writerow([f"resistance:{sp}", dg, c, t, round(c / t, 3) if t else 0])
        for gene, (d, t) in deletion.items():
            w.writerow(["deletion", gene, d, t, round(d / t, 3) if t else 0])
        for gene, (d, t) in cnv.items():
            w.writerow(["amplification", gene, d, t, round(d / t, 3) if t else 0])

    # HTML dashboard
    E = html.escape
    def bar(f):
        return f'<span style="display:inline-block;height:10px;width:{max(2,f*120):.0f}px;background:var(--blood);border-radius:5px;vertical-align:middle"></span>'
    H = [f"<!doctype html><html><head><meta charset=utf-8><title>{E(name)} surveillance</title>",
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600;700&display=swap">',
         "<style>:root{--paper:#F5F1EB;--ink:#231D19;--muted:#6C625A;--blood:#9E2B25;--teal:#0E6E66;--line:#DED4C5;--surface:#FCFAF6}"
         "body{background:var(--paper);color:var(--ink);font-family:'IBM Plex Sans',system-ui,sans-serif;padding:24px;line-height:1.5}"
         ".wrap{max-width:900px;margin:0 auto}h1{font:700 26px Georgia,serif}"
         "h2{font:700 12px monospace;letter-spacing:.06em;text-transform:uppercase;color:var(--blood);border-bottom:1px solid var(--line);padding-bottom:5px;margin:20px 0 8px}"
         "table{width:100%;border-collapse:collapse;font-size:12.5px}th{background:var(--ink);color:#fff;text-align:left;padding:5px 8px;font:600 10px monospace}"
         "td{padding:5px 8px;border-bottom:1px solid var(--line)}tr:nth-child(even) td{background:var(--surface)}"
         ".ey{font:600 11px monospace;letter-spacing:.16em;color:var(--blood);text-transform:uppercase}</style></head><body><div class=wrap>",
         "<div class=ey>MalaCall · molecular surveillance</div>",
         f"<h1>{E(name)}</h1><div style='color:var(--muted);font-size:13px'>{n} samples analysed</div>",
         "<h2>Species distribution</h2><table><tr><th>Species</th><th>Positive</th><th>Frequency</th><th></th></tr>"]
    for sp, c in sp_count.most_common():
        H.append(f"<tr><td>{E(sp)}</td><td>{c}/{n}</td><td>{c/n*100:.0f}%</td><td>{bar(c/n)}</td></tr>")
    H.append("</table><h2>Drug resistance frequency</h2><table><tr><th>Species</th><th>Drug</th><th>Resistant</th><th>Frequency</th><th></th></tr>")
    for sp in drug:
        for dg, c in sorted(drug[sp].items(), key=lambda x: -x[1]):
            t = drug_tot[sp][dg]
            if c:
                H.append(f"<tr><td>{E(sp)}</td><td>{E(dg)}</td><td>{c}/{t}</td><td>{c/t*100:.0f}%</td><td>{bar(c/t)}</td></tr>")
    H.append("</table><h2>RDT escape &amp; gene amplification</h2><table><tr><th>Marker</th><th>Type</th><th>Positive</th><th>Frequency</th></tr>")
    for gene, (d, t) in deletion.items():
        H.append(f"<tr><td>{E(gene)}</td><td>deletion (RDT escape)</td><td>{d}/{t}</td><td>{d/t*100:.0f}%</td></tr>")
    for gene, (d, t) in cnv.items():
        H.append(f"<tr><td>{E(gene)}</td><td>amplification</td><td>{d}/{t}</td><td>{d/t*100:.0f}%</td></tr>")
    H.append("</table><h2>Complexity of infection</h2><table><tr><th>Species</th><th>Samples</th><th>Mean COI</th><th>% polyclonal</th></tr>")
    for sp, vals in coi.items():
        poly = sum(1 for v in vals if v > 1)
        H.append(f"<tr><td>{E(sp)}</td><td>{len(vals)}</td><td>{sum(vals)/len(vals):.2f}</td><td>{poly/len(vals)*100:.0f}%</td></tr>")
    H.append("</table></div></body></html>")
    with open(f"{out}/{name}.cohort.html", "w") as fh:
        fh.write("\n".join(H))
    print(f"[cohort] {n} samples -> {out}/{name}.cohort.tsv + .html")


if __name__ == "__main__":
    main()
