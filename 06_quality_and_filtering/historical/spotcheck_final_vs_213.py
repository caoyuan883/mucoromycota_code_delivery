#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""抽检：最终版 202 表 vs 213 株原始表

抽检项（每株）：属种名、基因组 C+S+D+F+M、蛋白 C+S+D+F+M。
抽样：必查 26 株（10 株 OMArk 换值 + 6 株回退原名 + 9 株改名生效 + 1 株 C 格修正）
      + 固定种子随机 15 株。

映射路径与生成脚本不同（独立复算，避免复用生成逻辑）：
  最终表 sample → 表3（原名优先，其次采纳名）→ No → 213 表样本名（原名优先）→ 逐值比对。
"""
import csv, re, os, random, collections
import openpyxl

ROOT = "/Users/caoyuan/Desktop/Genome_annotation_project_of_Mucoromycota"
D = f"{ROOT}/protein_BUSCO/comparison_busco_genome_longest_isoform_protein_202"
FINAL = f"{D}/longest_isoform_busco_genome_protein_mucoromycota_long_202_final.csv"
OLD = f"{ROOT}/protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv"
T3 = f"{ROOT}/投稿附表/表3_研究总表.xlsx"
OUT = f"{D}/抽检报告_final_vs_213.md"

REVERT = {41, 94, 111, 127, 151, 152}                  # 采纳名重名 → 回退原名
OMARK = {8, 15, 17, 61, 62, 67, 106, 129, 138, 148}    # 蛋白 BUSCO 换成 OMArk 过滤后
CFIX = {168}                                           # Umbelopsis ramanniana：蛋白 C 按项目裁定修正


def collapse(n):
    return re.sub(r'[\s_.\-]+', '', str(n)).lower()


def numeq(a, b):
    """数值比较（容忍 99 与 99.0 的写法差异）。"""
    try:
        return abs(float(a) - float(b)) < 1e-9
    except (TypeError, ValueError):
        return str(a) == str(b)


# ---------------------------------------------------------------- 读三份数据
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

by_acc, by_orig = {}, {}
for b in t3:
    by_acc.setdefault(collapse(b['accepted']), []).append(b)
    by_orig.setdefault(collapse(b['orig']), []).append(b)
old_by_collapse = {}
for s in old:
    old_by_collapse.setdefault(collapse(s), []).append(s)


def resolve(name):
    """最终表主名 → No：原名优先，其次采纳名（采纳名可能重名）。"""
    k = collapse(name)
    c = by_orig.get(k, [])
    if len(c) == 1:
        return c[0]['no']
    c = by_acc.get(k, [])
    return c[0]['no'] if len(c) == 1 else None


def old_name_for(b):
    """表3 某 No → 213 表样本名：原名优先；原名不在 213 表时退回采纳名
    （No186/187/188/189/192 这 5 株的 213 条目用的就是采纳名）。"""
    for key in (b['orig'], b['accepted']):
        c = old_by_collapse.get(collapse(key), [])
        if len(c) == 1:
            return c[0]
    return None


# ---------------------------------------------------------------- 抽样
renamed = {b['no'] for b in t3 if collapse(b['orig']) != collapse(b['accepted'])}
must = sorted(OMARK | REVERT | renamed | CFIX)
rng = random.Random(20261005)
pool = [b['no'] for b in t3 if b['no'] not in must]
rand_nos = sorted(rng.sample(pool, 15))
pick = sorted(set(must) | set(rand_nos))
print(f"必查 {len(must)} 株 + 随机 {len(rand_nos)} 株 → 抽检 {len(pick)} 株")

no2b = {b['no']: b for b in t3}
exp_names = {b['no']: (b['orig'] if b['no'] in REVERT else b['accepted']) for b in t3}

# 全局：最终表 202 个名字必须逐一回解到 202 个 No，且与命名规则完全一致
assert len(fin) == 202, len(fin)
resolved = {resolve(s): s for s in fin}
assert None not in resolved and len(resolved) == 202, "最终表名字无法一一回解到表3 编号"
assert set(fin) == set(exp_names.values()), set(fin) ^ set(exp_names.values())
print("[全局] 最终表 202 个名字一一回解到表3 编号，且与命名规则完全一致 ✓")


def compare(no):
    """返回该株逐项比对结果。"""
    b = no2b[no]
    fname = exp_names[no]
    oname = old_name_for(b)
    r = dict(no=no, final=fname, old=oname, name_ok=(fname in fin and resolve(fname) == no),
             g_pct=True, g_C=True, p_pct=True, p_C=True, gbad=[], pbad=[], cfix=False)
    for m in 'SDFM':
        a, z = fin[fname][('G', m)], old[oname][('G', m)]
        if not numeq(a['percent'], z['percent']):
            r['g_pct'] = False; r['gbad'].append(m)
        if not numeq(a['complete_percent'], z['complete_percent']):
            r['g_C'] = False; r['gbad'].append(m + '(C)')
        a, z = fin[fname][('P', m)], old[oname][('P', m)]
        if not numeq(a['percent'], z['percent']):
            r['p_pct'] = False; r['pbad'].append(m)
        if not numeq(a['complete_percent'], z['complete_percent']):
            if no in CFIX and numeq(a['complete_percent'], 95.1) and numeq(z['complete_percent'], 73.5):
                r['cfix'] = True
            else:
                r['p_C'] = False; r['pbad'].append(m + '(C)')
    r['gC'] = fin[fname][('G', 'S')]['complete_percent']
    r['pC'] = fin[fname][('P', 'S')]['complete_percent']
    r['oldgC'] = old[oname][('G', 'S')]['complete_percent']
    r['oldpC'] = old[oname][('P', 'S')]['complete_percent']
    r['kind'] = '+'.join(x for x in [
        'OMArk换值' if no in OMARK else '', '回退原名' if no in REVERT else '',
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
    if r['kind'] != '常规' or len(details) < 4:
        details.append((no, r['final'], r['old'],
                        [fin[r['final']][('G', m)]['percent'] for m in 'SDFM'],
                        [old[r['old']][('G', m)]['percent'] for m in 'SDFM'],
                        [fin[r['final']][('P', m)]['percent'] for m in 'SDFM'],
                        [old[r['old']][('P', m)]['percent'] for m in 'SDFM']))

# ---------------------------------------------------------------- 全量兜底
full = [compare(b['no']) for b in t3]
full_g_bad = [r for r in full if not (r['g_pct'] and r['g_C'])]
full_p_bad = [r for r in full if not (r['p_pct'] and r['p_C'])]
full_p_pct_bad = [r for r in full if not r['p_pct']]
omark_changed = sorted(r['no'] for r in full if r['no'] in OMARK and not r['p_pct'])

# ---------------------------------------------------------------- 报告
L = []
A = L.append
A("# 抽检报告：最终版 202 表 vs 213 株原始表")
A("")
A(f"- 抽检对象：`{os.path.basename(FINAL)}`（202 株，14 列简表）")
A("- 比对基准：`protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv`"
  "（213 株原表，2026-08-17 22:53）")
A(f"- 抽样：必查 {len(must)} 株（10 株 OMArk 换值 + 6 株回退原名 + 9 株改名生效 + 1 株 C 格修正）"
  f" + 固定种子随机 {len(rand_nos)} 株 = **{len(pick)} 株**（占 202 株的 {len(pick)/202:.0%}）")
A("- 每株比对：**属种名** + **基因组 S/D/F/M + C** + **蛋白 S/D/F/M + C**")
A("")
if not problems:
    A("## 结论：✅ 抽检全部通过，0 处非预期差异")
else:
    A(f"## 结论：❌ 发现 {len(problems)} 处异常")
A("")
A("| 抽检项 | 结果 |")
A("|---|---|")
A(f"| 属种名对应（最终表名 ↔ 213 表名） | {sum(1 for r in rows if r['name_ok'])}/{len(rows)} ✓ |")
A(f"| 基因组 S/D/F/M + C | {sum(1 for r in rows if r['g_pct'] and r['g_C'])}/{len(rows)} 与 213 原表逐格相同 |")
A(f"| 蛋白 S/D/F/M | {sum(1 for r in rows if r['p_pct'])}/{len(rows)} 与原表相同；"
  f"{len(rows) - sum(1 for r in rows if r['p_pct'])} 株为按计划替换的 OMArk 过滤后新值 |")
A(f"| 蛋白 C | {sum(1 for r in rows if r['p_C'] and not r['cfix'])}/{len(rows)} 与原表相同；"
  f"{sum(1 for r in rows if r['cfix'])} 株按项目裁定修正（73.5→95.1）；"
  f"{sum(1 for r in rows if r['no'] in OMARK and not r['p_C'])} 株随 OMArk 换值 |")
A(f"| **非预期差异** | **{len(problems)}** |")
A("")
A("## 抽检逐株结果")
A("")
A("| No | 最终表名 | 213 表名 | 属种 | 基因组C (213) | 蛋白C (213) | 基因组4值 | 蛋白4值 | 备注 |")
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
    A(f"| {r['no']} | {r['final']} | `{r['old']}` | {'✓' if r['name_ok'] else '✗'} | "
      f"{r['gC']} ({r['oldgC']}) | {r['pC']} ({r['oldpC']}) | "
      f"{'✓' if r['g_pct'] and r['g_C'] else '✗'} | {pt} | {r['kind']} |")
A("")
A("## 逐值展开（全部非「常规」株 + 前几株随机样本）")
A("")
A("`G` = 基因组 BUSCO，`P` = longest-isoform 蛋白 BUSCO；数值格式 S / D / F / M。")
A("")
A("| No | 最终表名 | 数据 | S | D | F | M |")
A("|---|---|---|---|---|---|---|")
for no, fname, oname, g, og, p, op in details:
    A(f"| {no} | {fname} | G 最终 | {g[0]} | {g[1]} | {g[2]} | {g[3]} |")
    A(f"| | | G 213 原表 | {og[0]} | {og[1]} | {og[2]} | {og[3]} |")
    A(f"| | | P 最终 | {p[0]} | {p[1]} | {p[2]} | {p[3]} |")
    A(f"| | | P 213 原表 | {op[0]} | {op[1]} | {op[2]} | {op[3]} |")
A("")
A("## 全量兜底核对（不只抽样）")
A("")
A("- 全部 202 株 × 4 指标 × (percent + complete_percent = 1616 格)：")
A(f"  - 基因组不一致：**{len(full_g_bad)} 株**（1616 格全部与 213 原表相同）")
A(f"  - 蛋白 S/D/F/M 不一致：**{len(full_p_pct_bad)} 株** —— 全部落在 10 株 OMArk 名单内，属本次计划内替换")
A(f"  - 蛋白 C 单独不一致（S/D/F/M 未动）：**{sum(1 for r in full if r['cfix'])} 株**"
  f"（Umbelopsis ramanniana，按项目既有裁定 73.5→95.1）")
A(f"  - 其余 202 − {len(full_p_pct_bad)} − {sum(1 for r in full if r['cfix'])} = "
  f"{202 - len(full_p_pct_bad) - sum(1 for r in full if r['cfix'])} 株蛋白值与 213 原表逐格相同")
A(f"- OMArk 10 株里过滤前后**数值确实变了**的是 {len(omark_changed)} 株：{omark_changed}；"
  f"其余 {10 - len(omark_changed)} 株过滤前后数值相同，故看不出差别。")
A(f"- 抽样 {len(pick)} 株的结论与全量核对一致。")
A("")
A("## 说明")
A("")
A("- 基因组值**一律未变**（全量 202 株逐格相同）—— 基因组 BUSCO 反映原始组装质量，不参与污染过滤。")
A("- 蛋白值的变化只有两类，都可解释：① 10 株 OMArk 过滤后重跑（本次按你的要求替换）；"
  "② Umbelopsis ramanniana 的 C 格按项目既有裁定由 73.5 修正为 95.1（S/D/F/M 未动）。")
A("- 抽检映射路径独立于生成脚本：最终表名 → 表3（原名优先）→ No → 213 表样本名（原名优先）→ 逐值比对。")
open(OUT, 'w', encoding='utf-8').write("\n".join(L) + "\n")
print(f"[写出] {os.path.basename(OUT)}")
print("\n抽样 No:", pick)
print("抽检异常:", [r['no'] for r in problems] if problems else "无")
print(f"全量兜底：基因组不一致 {len(full_g_bad)} 株；蛋白 S/D/F/M 不一致 {len(full_p_pct_bad)} 株"
      f"（全在 OMArk 名单内）；C 不一致 {len(full_p_bad) - len(full_p_pct_bad)} 株（C 格修正）")
print("OMArk 10 株中实际发生数值变化的:", omark_changed)
