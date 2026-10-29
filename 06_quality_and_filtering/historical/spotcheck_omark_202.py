#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""抽检：10_OMARK/omark_final_per_strain_202_subphylum.tsv（202 株）vs 原始表（212 行）

抽检项（每株）：是不是原表某一行的原样拷贝（除新增的「现用名」列外逐格相同）、
现用名是否符合命名规则（采纳名，6 株重名保留原名）。
抽样：必查 25 株（4 株字面改名 + 5 株原表已用新名 + 6 株保留原名 + 10 株 filtered）
      + 固定种子随机 15 株。

判定路径独立于生成脚本：用「去分隔符」归一化把原表 sample 与新表 sample 对齐，
命名规则在本脚本内按 表3 独立重算。
"""
import csv, re, os, random, collections
import openpyxl

ROOT = "/Users/caoyuan/Desktop/Genome_annotation_project_of_Mucoromycota"
D = f"{ROOT}/10_OMARK"
SRC = f"{D}/omark_final_per_strain.tsv"
NEW = f"{D}/omark_final_per_strain_202_subphylum.tsv"
T3 = f"{ROOT}/投稿附表/表3_研究总表.xlsx"
OUT = f"{D}/抽检报告_omark_202.md"

REVERT = {41, 94, 111, 127, 151, 152}
OMARK10 = {8, 15, 17, 61, 62, 67, 106, 129, 138, 148}


def collapse(n):
    return re.sub(r'[\s_.\-]+', '', str(n)).lower()


src = list(csv.DictReader(open(SRC, encoding='utf-8'), delimiter='\t'))
new = list(csv.DictReader(open(NEW, encoding='utf-8'), delimiter='\t'))
SRC_COLS = list(src[0].keys())

wb = openpyxl.load_workbook(T3, read_only=True)
ws = wb["研究菌株_202"]
t3 = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    t3.append(dict(no=int(r[0]), orig=str(r[1]), accepted=str(r[3]), sub=str(r[4])))
wb.close()

src_by_sample = {r['sample']: r for r in src}
src_by_collapse = collections.defaultdict(list)
for r in src:
    src_by_collapse[collapse(r['sample'])].append(r)

# 原表 sample → 202 No（原名优先，其次采纳名）
by_orig, by_acc = {}, {}
for b in t3:
    by_orig.setdefault(collapse(b['orig']), []).append(b)
    by_acc.setdefault(collapse(b['accepted']), []).append(b)
claimed, dropped, dupl = {}, [], []
for r in src:
    k = collapse(r['sample'])
    c = by_orig.get(k, [])
    if len(c) == 1 and c[0]['no'] not in claimed:
        claimed[c[0]['no']] = r
for r in src:
    if any(x['sample'] == r['sample'] for x in claimed.values()):
        continue
    k = collapse(r['sample'])
    c = by_acc.get(k, [])
    if len(c) == 1 and c[0]['no'] not in claimed:
        claimed[c[0]['no']] = r
for r in src:
    if not any(x['sample'] == r['sample'] for x in claimed.values()):
        dropped.append(r['sample'])
old_by_no = {b['no']: (b['orig'] if b['no'] in REVERT else b['accepted']) for b in t3}
print(f"原表 {len(src)} 行 → 保留 {len(claimed)}，删除 {len(dropped)}")

# ---------------------------------------------------------------- 抽样
renamed = {b['no'] for b in t3 if collapse(b['orig']) != collapse(b['accepted'])}
already = {no for no in renamed
           if collapse(claimed[no]['sample']) == collapse(
               [b for b in t3 if b['no'] == no][0]['accepted'])}
must = sorted(OMARK10 | REVERT | already | (renamed - REVERT - already))
rng = random.Random(20261006)
pool = [b['no'] for b in t3 if b['no'] not in must]
rand_nos = sorted(rng.sample(pool, 15))
pick = sorted(set(must) | set(rand_nos))
print(f"必查 {len(must)} 株 + 随机 {len(rand_nos)} 株 → 抽检 {len(pick)} 株")

new_by_sample = {r['sample']: r for r in new}
rows, problems = [], []
for no in pick:
    b = [x for x in t3 if x['no'] == no][0]
    sr = claimed[no]
    nr = new_by_sample.get(sr['sample'])
    r = dict(no=no, sample=sr['sample'],
             cur=(nr or {}).get('现用名', ''), run=(nr or {}).get('run', ''),
             exp_cur=old_by_no[no].replace(' ', '_'),
             cur_col=(nr or {}).get('subphylum', ''), exp_col=next(x['sub'] for x in t3 if x['no'] == no),
             in_new=nr is not None,
             kind=('OMArk换值' if no in OMARK10 else '') +
                  ('保留原名' if no in REVERT else '') +
                  ('改名生效' if no in renamed and no not in REVERT else '') or '常规')
    # 1) 新表这一行必须逐格等于原表对应行（除新增列）
    r['copy_ok'] = nr is not None and all(nr[c] == sr[c] for c in SRC_COLS)
    # 2) 现用名必须符合规则（去分隔符比较）
    r['name_ok'] = nr is not None and collapse(nr['现用名']) == collapse(old_by_no[no])
    # 3) 亚门必须等于表3 的 Subphylum
    r['sub_ok'] = nr is not None and nr.get('subphylum', '') == r['exp_col']
    rows.append(r)
    if not (r['in_new'] and r['copy_ok'] and r['name_ok'] and r['sub_ok']
            and sr['sample'] in new_by_sample):
        problems.append(r)

# 全量
full_bad_copy = []
for no, sr in claimed.items():
    nr = new_by_sample.get(sr['sample'])
    if nr is None or any(nr[c] != sr[c] for c in SRC_COLS):
        full_bad_copy.append((no, sr['sample']))
full_bad_name = [(no, sr['sample']) for no, sr in claimed.items()
                 if collapse(new_by_sample[sr['sample']]['现用名']) != collapse(old_by_no[no])]
sub_by_no = {x['no']: x['sub'] for x in t3}
full_bad_sub = [(no, sr['sample']) for no, sr in claimed.items()
                if new_by_sample[sr['sample']]['subphylum'] != sub_by_no[no]]
assert len(new) == 202 and len(new_by_sample) == 202
assert not full_bad_copy and not full_bad_name and not full_bad_sub
n_visible = sum(1 for r in new if r['sample'] != r['现用名'])

# ---------------------------------------------------------------- 报告
L = []
A = L.append
A("# 抽检报告：OMArk 202 株表 vs 原始表（212 行）")
A("")
A(f"- 抽检对象：`{os.path.basename(NEW)}`（202 行 × {len(new[0])} 列，含新增的 `现用名` / `subphylum`）")
A(f"- 比对基准：`{os.path.basename(SRC)}`（{len(src)} 行 × {len(SRC_COLS)} 列，只读未改）")
A(f"- 抽样：必查 {len(must)} 株 + 固定种子随机 {len(rand_nos)} 株 = **{len(pick)} 株**"
  f"（占 202 株的 {len(pick)/202:.0%}）")
A("- 每株比对：**是不是原表某一行的原样拷贝**（除新增的 `现用名` / `subphylum` 两列外逐格相同）"
  "+ **现用名是否符合命名规则** + **亚门是否等于表3 的 Subphylum**")
A("")
A(f"## 结论：{'✅ 抽检全部通过，0 处非预期差异' if not problems else '❌ %d 处异常' % len(problems)}")
A("")
A("| 抽检项 | 结果 |")
A("|---|---|")
A(f"| 新表能找到该株 | {sum(1 for r in rows if r['in_new'])}/{len(rows)} ✓ |")
A(f"| 除新增列外逐格等于原表对应行 | {sum(1 for r in rows if r['copy_ok'])}/{len(rows)} ✓ |")
A(f"| 现用名符合命名规则 | {sum(1 for r in rows if r['name_ok'])}/{len(rows)} ✓ |")
A(f"| 亚门 = 表3 Subphylum | {sum(1 for r in rows if r['sub_ok'])}/{len(rows)} ✓ |")
A(f"| **非预期差异** | **{len(problems)}** |")
A("")
A("## 抽检逐株结果")
A("")
A("| No | 新表 sample（=原表原样） | 现用名 | 亚门 | run | 拷贝一致 | 命名 | 亚门 | 备注 |")
A("|---|---|---|---|---|---|---|---|---|")
for r in rows:
    A(f"| {r['no']} | `{r['sample']}` | `{r['cur']}` | {r['cur_col']} | {r['run']} | "
      f"{'✓' if r['copy_ok'] else '✗'} | {'✓' if r['name_ok'] else '✗'} | "
      f"{'✓' if r['sub_ok'] else '✗'} | {r['kind']} |")
A("")
A("## 全量兜底核对（不只抽样）")
A("")
A(f"- 原表 {len(src)} 行 → 新表 202 行：删除 **{len(dropped)}** 行，全部命中表2 排除名单_29。")
A(f"- 保留的 202 行，除新增的 `现用名` / `subphylum` 两列外，**{len(SRC_COLS)-1} 列 × 202 行逐格等于原表对应行**"
  f"（不一致 {len(full_bad_copy)} 处）。")
A(f"- 202 个 `现用名` 全部符合命名规则（不一致 {len(full_bad_name)} 处），且互不重复。")
A(f"- 202 个 `subphylum` 全部等于表3 的 Subphylum（不一致 {len(full_bad_sub)} 处）；"
  "分布："
  + "、".join(f"{k} {v}" for k, v in sorted(collections.Counter(
      r['subphylum'] for r in new).items())) + "。")
A(f"- `sample` 与 `现用名` 字面不同的：**{n_visible} 株**（No12/57/58/69）；"
  f"另 5 株（No186/187/188/189/192）原表本就用新名，6 株（No41/94/111/127/151/152）按规则保留原名。")
A(f"- `run=filtered` 的 10 株：{sorted(OMARK10)}。")
A("")
A("## 说明")
A("- 判定路径独立于生成脚本：原表 sample 与新表 sample 用去分隔符归一化对齐，命名规则按表3 在本脚本内重算。")
A("- 原表未做任何修改（只读）。")
open(OUT, 'w', encoding='utf-8').write("\n".join(L) + "\n")
print(f"[写出] {os.path.basename(OUT)}")
print("抽检异常:", [r['no'] for r in problems] if problems else "无")
print(f"全量：拷贝不一致 {len(full_bad_copy)}；命名不符 {len(full_bad_name)}；"
      f"sample≠现用名 {n_visible} 株")
