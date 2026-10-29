#!/usr/bin/env python3
import sys, re, difflib
from pathlib import Path
from collections import Counter

LIST   = Path(sys.argv[1]) if len(sys.argv)>1 else None
ING_DIR = Path("/home/py/Desktop/Mucoromycota_protein_fungi_odb12_single_copy_BUSCO")
ING_IDS = Path("/home/py/Desktop/5_outgroup/ingroup_protein_locus_ids.tsv")
OUT_IDS = Path("/home/py/Desktop/5_outgroup/single_copy_locus_ids.tsv")
REPORT  = Path("/home/py/Desktop/5_outgroup/list_verification.tsv")

if not LIST or not LIST.is_file():
    sys.exit("用法: verify_species_list.py <202株名单.txt>   （每行一个物种名）")

def norm(s):  # 仅用于比对，报告仍用原名
    return re.sub(r"[^a-z0-9]", "", s.lower())

wanted=[l.strip() for l in LIST.read_text(encoding="utf-8-sig").splitlines() if l.strip()]
dupes=[k for k,v in Counter(map(norm,wanted)).items() if v>1]
print(f"名单条目：{len(wanted)}   规范化后唯一：{len(set(map(norm,wanted)))}")
if dupes: print(f"  ⚠ 名单内有重复：{dupes}")

have_dir = {p.name for p in ING_DIR.iterdir() if p.is_dir()} if ING_DIR.is_dir() else set()
per={}
for line in ING_IDS.read_text().splitlines()[1:]:
    f=line.split("\t")
    if len(f)>=4: per[f[1]]=set(f[3].split(","))
print(f"已提取目录：{len(have_dir)}   位点清单条目：{len(per)}")
if have_dir != set(per):
    print(f"  ⚠ 目录与清单不一致：仅目录 {sorted(have_dir-set(per))[:5]}  仅清单 {sorted(set(per)-have_dir)[:5]}")

nmap={norm(s):s for s in per}
matched={}; missing=[]
for w in wanted:
    k=norm(w)
    if k in nmap: matched[w]=nmap[k]
    else: missing.append(w)
extra=[s for s in per if norm(s) not in {norm(w) for w in wanted}]

print(f"\n=== 核对结果 ===")
print(f"  精确匹配（含大小写/标点差异）：{len(matched)}")
print(f"  名单中有、结果里找不到：{len(missing)}")
print(f"  结果里有、名单外：{len(extra)}")

if missing:
    print("\n  名单中缺失的条目及最相近候选：")
    pool=list(per)
    for w in missing:
        cand=difflib.get_close_matches(w, pool, n=3, cutoff=0.75)
        print(f"    {w:38s} -> {cand if cand else '无相近项'}")
if extra:
    print(f"\n  结果里多出（应为被剔除的株）共 {len(extra)}：")
    for s in sorted(extra): print(f"    {s}")

with REPORT.open("w") as fh:
    fh.write("list_name\tmatched_dir_name\tstatus\tn_loci\n")
    for w in wanted:
        if w in matched:
            d=matched[w]; fh.write(f"{w}\t{d}\tMATCH\t{len(per[d])}\n")
        else:
            fh.write(f"{w}\t\tMISSING\t\n")
    for s in sorted(extra):
        fh.write(f"\t{s}\tEXTRA_NOT_IN_LIST\t{len(per[s])}\n")
print(f"\n核对表 -> {REPORT}")

if missing:
    print("\n⚠ 存在缺失，占有率统计暂不输出。请先修正名单或目录名后重跑。")
    sys.exit(0)

sub={matched[w]:per[matched[w]] for w in wanted}
N=len(sub); occ=Counter()
for ids in sub.values():
    for b in ids: occ[b]+=1
print(f"\n=== 内群 {N} 株 占有率 -> 位点数 ===")
for t in (1.00,0.99,0.98,0.95,0.90,0.85,0.80):
    print(f"   >= {t*100:5.1f}%  ({int(t*N):3d}/{N})   {sum(1 for c in occ.values() if c>=t*N):5d}")

og={}
if OUT_IDS.is_file():
    for line in OUT_IDS.read_text().splitlines()[1:]:
        f=line.split("\t")
        if len(f)>=4 and f[0]=="protein": og[f[1]]=set(f[3].split(","))
if og:
    print(f"\n外类群（蛋白模式）{len(og)} 株：", ", ".join(f"{k}={len(v)}" for k,v in sorted(og.items())))
    allsp={**sub,**og}; M=len(allsp); occ2=Counter()
    for ids in allsp.values():
        for b in ids: occ2[b]+=1
    og_all=set.intersection(*og.values())
    print(f"\n=== 全 {M} 末端 占有率 -> 位点数（蛋白模式）===")
    for t in (1.00,0.99,0.98,0.95,0.90,0.85,0.80):
        n_all=sum(1 for c in occ2.values() if c>=t*M)
        n_og =sum(1 for b,c in occ2.items() if c>=t*M and b in og_all)
        print(f"   >= {t*100:5.1f}%  ({int(t*M):3d}/{M})   共 {n_all:5d}   其中5外类群全有 {n_og:5d}")
    print("\n  ★ 建议以'5 外类群全有'那一列为准——定根质量取决于外类群的数据完整性。")
else:
    print("\n（未找到外类群蛋白模式位点清单，跳过 207 末端统计）")
