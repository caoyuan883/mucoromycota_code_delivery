#!/usr/bin/env python3
import re
from pathlib import Path
ROOT = Path.home()/"bioinfo/benchmark12_busco_ref_protein/results_mucoromycota"
OUT  = Path("/home/py/Desktop/12_faa_12_gff_technical_validation/busco_reference_12.tsv")
pat = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")
rows=[]
for d in sorted(ROOT.iterdir()):
    if not d.is_dir(): continue
    ss = next(iter(sorted(d.glob("short_summary*.txt"))), None)
    if ss is None: print(f"[缺summary] {d.name}"); continue
    m = pat.search(ss.read_text())
    if not m: print(f"[解析失败] {d.name}"); continue
    C,S,D,F,M,n = m.groups()
    rows.append((d.name,C,S,D,F,M,n))
    print(f"{d.name:32s} C={C:6s} S={S:6s} D={D:5s} F={F:5s} M={M:5s} n={n}")
with OUT.open("w") as fh:
    fh.write("sample\tC\tS\tD\tF\tM\tn\n")
    for r in rows: fh.write("\t".join(r)+"\n")
print(f"\n{len(rows)} 株 -> {OUT}")
