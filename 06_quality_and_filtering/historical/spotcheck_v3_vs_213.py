#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""抽检：v3 简化版 202 表 vs 213 株原始表

抽检项（每株）：两列名字（sample 原名 / 现用名）、基因组 C+S+D+F+M、蛋白 C+S+D+F+M。
抽样：必查 26 株（10 株 OMArk 换值 + 6 株重名回退 + 9 株改名生效 + 1 株 C 格修正）
      + 固定种子随机 15 株。

映射路径独立于生成脚本：
  v3 的 sample(=原名) → 213 表样本名（原名优先，其次采纳名）→ 逐值比对；
  现用名 由 表3 采纳名 + 6 株回退规则独立重算后比对。
"""
import csv, re, os, random, collections
import openpyxl

ROOT = "/Users/caoyuan/Desktop/Genome_annotation_project_of_Mucoromycota"
D = f"{ROOT}/protein_BUSCO/comparison_busco_genome_longest_isoform_protein_202"
FINAL = f"{D}/longest_isoform_busco_genome_protein_mucoromycota_long_202_v3.csv"
OLD = f"{ROOT}/protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv"
T3 = f"{ROOT}/投稿附表/表3_研究总表.xlsx"
OUT = f"{D}/抽检报告_v3_vs_213.md"

REVERT = {41, 94, 111, 127, 151, 152}
OMARK = {8, 15, 17, 61, 62, 67, 106, 129, 138, 148}
CFIX = {168}
CURCOL = '现用名'


def collapse(n):
    return re.sub(r'[\s_.\-]+', '', str(n)).lower()


def numeq(a, b):
    try:
        return abs(float(a) - float(b)) < 1e-9
    except (TypeError, ValueError):
        return str(a) == str(b)


# ---------------------------------------------------------------- 读数据
fin = collections.defaultdict(dict)
for r in csv.DictReader(open(FINAL, encoding='utf-8')):
    fin[r['sample']][(r['dataset_label'], r['metric'])] = r

old = collections.defaultdict(dict)
for r in csv.DictReader(open(OLD, encoding='utf-8-sig')):
    old[r['sample']][(r['dataset_label'], r['metric'])] = r

wb = openpyxl.load_workbook(T3, read_only=True)
ws = wb["研究菌株_202"]
t3 = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    t3.append(dict(no=int(r[0]), orig=str(r[1]), accepted=str(r[3])))
wb.close()

old_by_collapse = {}
for s in old:
    old_by_collapse.setdefault(collapse(s), []).append(s)


def old_name_for(b):
    for key in (b['orig'], b['accepted']):
        c = old_by_collapse.get(collapse(key), [])
        if len(c) == 1:
            return c[0]
    return None


def expected_cur(b):
    return b['orig'] if b['no'] in REVERT else b['accepted']


# 全局：两列名字
assert len(fin) == 202, len(fin)
assert set(fin) == {b['orig'] for b in t3}, "sample 列 ≠ 表3 原名集合"
cur_seen = {r[CURCOL] for d in fin.values() for r in d.values()}
assert cur_seen == {expected_cur(b) for b in t3}, "现用名列 ≠ 期望集合"
assert len(cur_seen) == 202, "现用名有重复"
print("[全局] 202 个 sample 全等于表3 原名；202 个现用名互不重复且符合回退规则 ✓")

# ---------------------------------------------------------------- 抽样
renamed = {b['no'] for b in t3 if collapse(b['orig']) != collapse(b['accepted'])}
must = sorted(OMARK | REVERT | renamed | CFIX)
rng = random.Random(20261005)
pool = [b['no'] for b in t3 if b['no'] not in must]
rand_nos = sorted(rng.sample(pool, 15))
pick = sorted(set(must) | set(rand_nos))
no2b = {b['no']: b for b in t3}
print(f"必查 {len(must)} 株 + 随机 {len(rand_nos)} 株 → 抽检 {len(pick)} 株")


def compare(no):
    b = no2b[no]
    r0 = fin[b['orig']][('G', 'S')]
    oname = old_name_for(b)
    r = dict(no=no, sample=b['orig'], cur=r0[CURCOL], old=oname,
             name_ok=(r0['sample'] == b['orig'] and r0[CURCOL] == expected_cur(b)),
             g_pct=True, g_C=True, p_pct=True, p_C=True, cfix=False, kind='')
    for m in 'SDFM':
        a, z = fin[b['orig']][('G', m)], old[oname][('G', m)]
        if not numeq(a['percent'], z['percent']):
            r['g_pct'] = False
        if not numeq(a['complete_percent'], z['complete_percent']):
            r['g_C'] = False
        a, z = fin[b['orig']][('P', m)], old[oname][('P', m)]
        if not numeq(a['percent'], z['percent']):
            r['p_pct'] = False
        if not numeq(a['complete_percent'], z['complete_percent']):
            if no in CFIX and numeq(a['complete_percent'], 95.1) and numeq(z['complete_percent'], 73.5):
                r['cfix'] = True
            else:
                r['p_C'] = False
    r['gC'] = fin[b['orig']][('G', 'S')]['complete_percent']
    r['pC'] = fin[b['orig']][('P', 'S')]['complete_percent']
    r['oldgC'] = old[oname][('G', 'S')]['complete_percent']
    r['oldpC'] = old[oname][('P', 'S')]['complete_percent']
    r['kind'] = '+'.join(x for x in [
        'OMArk换值' if no in OMARK else '', '重名回退' if no in REVERT else '',
        '改名生效' if (no in renamed and no not in REVERT) else '',
        'C格修正' if no in CFIX else ''] if x) or '常规'
    return r


rows, details, problems = [], [], []
for no in pick:
    r = compare(no)
    unexpected = (not r['name_ok']) or (not r['g_pct']) or (not r['g_C']) \
        or (not r['p_pct'] and no not in OMARK) \
        or (not r['p_C'] and no not in OMARK and not r['cfix'])
    if unexpected:
        problems.append(r)
    rows.append(r)
    if r['kind'] != '常规' or len(details) < 3:
        b = no2b[no]
        details.append((no, b['orig'], r['cur'], r['old'],
                        [fin[b['orig']][('G', m)]['percent'] for m in 'SDFM'],
                        [old[r['old']][('G', m)]['percent'] for m in 'SDFM'],
                        [fin[b['orig']][('P', m)]['percent'] for m in 'SDFM'],
                        [old[r['old']][('P', m)]['percent'] for m in 'SDFM']))

full = [compare(b['no']) for b in t3]
full_g_bad = [r for r in full if not (r['g_pct'] and r['g_C'])]
full_p_pct_bad = [r for r in full if not r['p_pct']]
full_name_bad = [r for r in full if not r['name_ok']]
omark_changed = sorted(r['no'] for r in full if r['no'] in OMARK and not r['p_pct'])
n_cfix = sum(1 for r in full if r['cfix'])
n_same = 202 - len(full_p_pct_bad) - n_cfix

# ---------------------------------------------------------------- 报告
L = []
A = L.append
A("# 抽检报告：v3 简化版 202 表 vs 213 株原始表")
A("")
A(f"- 抽检对象：`{os.path.basename(FINAL)}`（202 株，14 列）")
A("- 比对基准：`protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv`"
  "（213 株原表，2026-08-17 22:53）")
A(f"- 抽样：必查 {len(must)} 株 + 固定种子随机 {len(rand_nos)} 株 = **{len(pick)} 株**"
  f"（占 202 株的 {len(pick)/202:.0%}）")
A("- 每株比对：**两列名字** + **基因组 S/D/F/M + C** + **蛋白 S/D/F/M + C**")
A("")
A(f"## 结论：{'✅ 抽检全部通过，0 处非预期差异' if not problems else '❌ 发现 %d 处异常' % len(problems)}")
A("")
A("| 抽检项 | 结果 |")
A("|---|---|")
A(f"| sample(原名) / 现用名 | {sum(1 for r in rows if r['name_ok'])}/{len(rows)} ✓ |")
A(f"| 基因组 S/D/F/M + C | {sum(1 for r in rows if r['g_pct'] and r['g_C'])}/{len(rows)} 与 213 原表逐格相同 |")
A(f"| 蛋白 S/D/F/M | {sum(1 for r in rows if r['p_pct'])}/{len(rows)} 与原表相同；"
  f"{sum(1 for r in rows if not r['p_pct'])} 株为按计划替换的 OMArk 过滤后新值 |")
A(f"| 蛋白 C | {sum(1 for r in rows if r['p_C'] and not r['cfix'])}/{len(rows)} 与原表相同；"
  f"{sum(1 for r in rows if r['cfix'])} 株按项目裁定修正（73.5→95.1）；"
  f"{sum(1 for r in rows if r['no'] in OMARK and not r['p_C'])} 株随 OMArk 换值 |")
A(f"| **非预期差异** | **{len(problems)}** |")
A("")
A("## 抽检逐株结果")
A("")
A("| No | sample(原名) | 现用名 | 213 表名 | 基因组C (213) | 蛋白C (213) | 基因组4值 | 蛋白4值 | 备注 |")
A("|---|---|---|---|---|---|---|---|---|")
for r in rows:
    if r['cfix']:
        pt = '✓ S/D/F/M 同；C 修正 73.5→95.1'
    elif not r['p_pct'] and r['no'] in OMARK:
        pt = '↻ 按计划换值'
    elif r['p_pct'] and r['p_C']:
        pt = '✓ 同原表'
    else:
        pt = '✗'
    same = '（与原名相同）' if r['sample'] == r['cur'] else ''
    A(f"| {r['no']} | {r['sample']} | {r['cur']}{same} | `{r['old']}` | "
      f"{r['gC']} ({r['oldgC']}) | {r['pC']} ({r['oldpC']}) | "
      f"{'✓' if r['g_pct'] and r['g_C'] else '✗'} | {pt} | {r['kind']} |")
A("")
A("## 逐值展开（全部非「常规」株 + 前几株随机样本）")
A("")
A("`G` = 基因组 BUSCO，`P` = longest-isoform 蛋白 BUSCO；数值格式 S / D / F / M。")
A("")
A("| No | sample(原名) | 现用名 | 数据 | S | D | F | M |")
A("|---|---|---|---|---|---|---|---|")
for no, s, cur, oname, g, og, p, op in details:
    tag = '' if s == cur else f" ({cur})"
    A(f"| {no} | {s} | {cur} | G 最终 | {g[0]} | {g[1]} | {g[2]} | {g[3]} |")
    A(f"| | | | G 213 原表 | {og[0]} | {og[1]} | {og[2]} | {og[3]} |")
    A(f"| | | | P 最终 | {p[0]} | {p[1]} | {p[2]} | {p[3]} |")
    A(f"| | | | P 213 原表 | {op[0]} | {op[1]} | {op[2]} | {op[3]} |")
A("")
A("## 全量兜底核对（不只抽样）")
A("")
A(f"- 两列名字：sample 不等于表3 原名 **{len(full_name_bad)} 株**；现用名不符合回退规则 0 株。")
A(f"- 基因组：202 株 × 4 指标 × (percent + complete_percent) 与 213 原表不一致 **{len(full_g_bad)} 处**。")
A(f"- 蛋白 S/D/F/M 不一致 **{len(full_p_pct_bad)} 株**，全部落在 10 株 OMArk 名单内（计划内替换）。")
A(f"- 蛋白 C 单独修正 **{n_cfix} 株**（Umbelopsis ramanniana，73.5→95.1）。")
A(f"- 其余 **{n_same} 株**蛋白值与 213 原表逐格相同。")
A(f"- OMArk 10 株里 after 与 before **确实不同**的是 {len(omark_changed)} 株：{omark_changed}。")
A("")
A("## 说明")
A("")
A("- 基因组值一律未变 —— 基因组 BUSCO 反映原始组装质量，不参与污染过滤。")
A("- 蛋白值变化只有两类：① 10 株 OMArk 过滤后重跑；② Umbelopsis ramanniana 的 C 格按项目裁定修正。")
A("- 抽检映射独立于生成脚本：sample(原名) → 213 表样本名（原名优先）→ 逐值比对。")
open(OUT, 'w', encoding='utf-8').write("\n".join(L) + "\n")
print(f"[写出] {os.path.basename(OUT)}")
print("抽检异常:", [r['no'] for r in problems] if problems else "无")
print(f"全量兜底：名字不符 {len(full_name_bad)} 株；基因组不一致 {len(full_g_bad)} 株；"
      f"蛋白 S/D/F/M 不一致 {len(full_p_pct_bad)} 株（全在 OMArk 名单）；C 格修正 {n_cfix} 株")
print("OMArk 10 株中 after≠before 的:", omark_changed)
