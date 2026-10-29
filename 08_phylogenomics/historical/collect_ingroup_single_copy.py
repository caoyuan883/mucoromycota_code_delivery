#!/usr/bin/env python3
import shutil, sys
from pathlib import Path

SRC  = Path("/home/py/bioinfo/mucoromycota_busco_longest_isoform_213/results_fungi")
DST  = Path("/home/py/Desktop/Mucoromycota_protein_fungi_odb12_single_copy_BUSCO")
IDS  = Path("/home/py/Desktop/5_outgroup/ingroup_protein_locus_ids.tsv")
SUF  = "_longest_isoform"
WHITE = Path(sys.argv[1]) if len(sys.argv) > 1 else None   # 可选：202株名单，每行一个物种名

keep = None
if WHITE and WHITE.is_file():
    keep = {l.strip() for l in WHITE.read_text().splitlines() if l.strip()}
    print(f"启用白名单：{len(keep)} 株\n")

DST.mkdir(parents=True, exist_ok=True)
per={}; skipped=[]; issues=[]
dirs = sorted(p for p in SRC.iterdir() if p.is_dir())
print(f"源目录物种数：{len(dirs)}\n")

for i, sp_dir in enumerate(dirs, 1):
    name = sp_dir.name
    sample = name[:-len(SUF)] if name.endswith(SUF) else name
    if keep is not None and sample not in keep:
        skipped.append(sample); continue
    sc = next(iter(sp_dir.rglob("single_copy_busco_sequences")), None)
    if sc is None:
        print(f"  [缺] {sample}: 无 single_copy_busco_sequences"); continue
    out = DST/sample; out.mkdir(exist_ok=True)
    ids=set(); n_multi=0; n_star=0
    for faa in sc.glob("*.faa"):
        shutil.copy2(faa, out/faa.name)
        ids.add(faa.stem)
        t = faa.read_text()
        if t.count(">")>1: n_multi+=1
        if "*" in t: n_star+=1
    per[sample]=ids
    if n_multi or n_star: issues.append((sample,n_multi,n_star))
    if i % 20 == 0 or i == len(dirs):
        print(f"  [{i}/{len(dirs)}] 已处理 {len(per)} 株 …")

print(f"\n完成：{len(per)} 株 -> {DST}")
if skipped: print(f"白名单外跳过：{len(skipped)} 株")

cnt=sorted((len(v),k) for k,v in per.items())
print(f"\n位点数  中位 {cnt[len(cnt)//2][0]}   最少 {cnt[0][0]} ({cnt[0][1]})   最多 {cnt[-1][0]} ({cnt[-1][1]})")
print("最少的 8 株：")
for n,s in cnt[:8]: print(f"   {s:34s} {n}")

inter=set.intersection(*per.values()); union=set.union(*per.values())
print(f"\n全部 {len(per)} 株共有位点：{len(inter)}    并集：{len(union)}")
from collections import Counter
occ=Counter()
for ids in per.values():
    for b in ids: occ[b]+=1
N=len(per)
print("\n占有率阈值 -> 位点数（仅内群）")
for th in (1.00,0.95,0.90,0.85,0.80,0.75):
    print(f"   >= {th*100:5.1f}%  ({int(th*N):3d}/{N} 株)   {sum(1 for c in occ.values() if c>=th*N):5d}")

with IDS.open("w") as fh:
    fh.write("mode\tsample\tn_loci\tbusco_ids\n")
    for s,ids in sorted(per.items()):
        fh.write(f"protein\t{s}\t{len(ids)}\t{','.join(sorted(ids))}\n")
print(f"\n位点清单 -> {IDS}")

if issues:
    print(f"\n需留意（多序列文件 / 含*号）共 {len(issues)} 株，前 10：")
    for s,m,st in issues[:10]: print(f"   {s:34s} 多序列 {m:5d}   含* {st:5d}")
else:
    print("\n全部文件均为单序列且不含 * ，与外类群蛋白模式一致 ✓")
