#!/usr/bin/env python3
from pathlib import Path
from collections import Counter
import sys

SRC = Path("/home/py/Desktop/5_outgroup/ingroup_protein_locus_ids.tsv")
OUT = Path("/home/py/Desktop/5_outgroup/species_list_213_protein_loci.tsv")
WHITE = Path(sys.argv[1]) if len(sys.argv) > 1 else None

rows=[]; per={}
for line in SRC.read_text().splitlines()[1:]:
    f=line.split("\t")
    if len(f)<4: continue
    per[f[1]]=set(f[3].split(","))
    rows.append((f[1], int(f[2])))

keep=None
if WHITE and WHITE.is_file():
    keep={l.strip() for l in WHITE.read_text().splitlines() if l.strip()}

rows.sort(key=lambda x: x[1])
with OUT.open("w") as fh:
    fh.write("rank\tsample\tn_loci\tin_202\n")
    for i,(s,n) in enumerate(sorted(rows,key=lambda x:x[0]),1):
        flag = ("Y" if s in keep else "N") if keep else "?"
        fh.write(f"{i}\t{s}\t{n}\t{flag}\n")
print(f"{len(rows)} 株 -> {OUT}\n")

print("位点数最少的 15 株：")
for s,n in rows[:15]: print(f"   {s:36s} {n}")

if keep:
    sub={k:v for k,v in per.items() if k in keep}
    missing = keep - set(per)
    extra   = set(per) - keep
    print(f"\n名单比对：匹配 {len(sub)}  名单中缺失 {len(missing)}  多余 {len(extra)}")
    if missing: print("  名单中有但结果里没有：", sorted(missing)[:20])
    if extra:   print("  结果里有但名单外：", sorted(extra)[:20])
    if sub:
        N=len(sub); occ=Counter()
        for ids in sub.values():
            for b in ids: occ[b]+=1
        print(f"\n占有率阈值 -> 位点数（{N} 株口径）")
        for t in (1.00,0.98,0.95,0.90,0.85,0.80):
            print(f"   >= {t*100:5.1f}%  ({int(t*N):3d}/{N})   {sum(1 for c in occ.values() if c>=t*N):5d}")
