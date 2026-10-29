#!/usr/bin/env python3
import re
from pathlib import Path
BASE=Path("/home/py/Desktop/5_outgroup")
pat=re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")
rows=[]
for mode,sub in [("genome","busco_genome"),("protein","busco_protein")]:
    root=BASE/sub
    if not root.is_dir(): continue
    for d in sorted(root.iterdir()):
        if not d.is_dir(): continue
        ss=next(iter(sorted(d.glob("short_summary*.txt"))),None)
        if ss is None: print(f"[缺summary] {mode} {d.name}"); continue
        m=pat.search(ss.read_text())
        if not m: print(f"[解析失败] {mode} {d.name}"); continue
        sc=next(iter(d.rglob("single_copy_busco_sequences")),None)
        n=len(list(sc.glob("*.faa"))) if sc else 0
        rows.append((d.name,mode,*m.groups(),n))
        print(f"{d.name:26s} {mode:8s} C={m.group(1):6s} S={m.group(2):6s} "
              f"D={m.group(3):5s} F={m.group(4):5s} M={m.group(5):6s} n={m.group(6)}  单拷贝faa={n}")
out=BASE/"outgroup_busco_summary.tsv"
with out.open("w") as fh:
    fh.write("sample\tmode\tC\tS\tD\tF\tM\tn\tn_single_copy_faa\n")
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
print(f"\n{len(rows)} 条 -> {out}")
