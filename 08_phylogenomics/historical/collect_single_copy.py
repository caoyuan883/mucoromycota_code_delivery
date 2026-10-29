#!/usr/bin/env python3
import shutil
from pathlib import Path

BASE = Path("/home/py/Desktop/5_outgroup")
JOBS = [("genome",  BASE/"busco_genome",  BASE/"genome_single_copy"),
        ("protein", BASE/"busco_protein", BASE/"protein_single_copy")]

report={}; oddities=[]
for mode, src, dst in JOBS:
    if not src.is_dir():
        print(f"[跳过] {mode}: 未找到 {src}"); continue
    dst.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {mode} -> {dst} ===")
    per={}
    for sp_dir in sorted(p for p in src.iterdir() if p.is_dir()):
        sp = sp_dir.name
        sc = next(iter(sp_dir.rglob("single_copy_busco_sequences")), None)
        if sc is None:
            print(f"  [缺] {sp}"); continue
        out = dst/sp; out.mkdir(exist_ok=True)
        ids=set(); n_multi=0; n_star=0
        for faa in sorted(sc.glob("*.faa")):
            shutil.copy2(faa, out/faa.name)          # 原样复制，不动序列头
            ids.add(faa.stem)
            txt=faa.read_text()
            if txt.count(">")>1: n_multi+=1
            if "*" in txt: n_star+=1
        per[sp]=ids
        print(f"  {sp:26s} 位点 {len(ids):5d}   多序列文件 {n_multi:3d}   含*号 {n_star:5d}")
        if n_multi: oddities.append((mode,sp,"多序列文件",n_multi))
    report[mode]=per
    if per:
        print(f"  5株共有位点: {len(set.intersection(*per.values()))}"
              f"   并集: {len(set.union(*per.values()))}")

if "genome" in report and "protein" in report:
    print("\n=== 同一物种 基因组模式 vs 蛋白模式 位点重合 ===")
    for sp in sorted(set(report["genome"]) & set(report["protein"])):
        g,p = report["genome"][sp], report["protein"][sp]
        print(f"  {sp:26s} genome {len(g):5d}  protein {len(p):5d}  "
              f"交集 {len(g&p):5d}  仅genome {len(g-p):4d}  仅protein {len(p-g):4d}")
    gi=set.intersection(*report["genome"].values())
    pi=set.intersection(*report["protein"].values())
    print(f"\n  ★ 5株 × 两模式 全部具备的位点: {len(gi & pi)}")

with (BASE/"single_copy_locus_ids.tsv").open("w") as fh:
    fh.write("mode\tsample\tn_loci\tbusco_ids\n")
    for mode,d in report.items():
        for sp,ids in d.items():
            fh.write(f"{mode}\t{sp}\t{len(ids)}\t{','.join(sorted(ids))}\n")
print(f"\n位点清单 -> {BASE/'single_copy_locus_ids.tsv'}")
if oddities:
    print("\n需留意:"); [print("  ",*o) for o in oddities]
