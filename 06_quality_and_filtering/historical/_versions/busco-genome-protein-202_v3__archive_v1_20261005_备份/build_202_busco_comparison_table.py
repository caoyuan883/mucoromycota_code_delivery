#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 213 株的 BUSCO 基因组–蛋白对比表（mucoromycota_odb12）重建为当前审计的 202 株版本。

原则
----
* 只读。旧表（213）一个字节都不改；本脚本只读它并把行「拷贝」到新表。
* 主键 = 表3_研究总表.xlsx「研究菌株_202」的权威编号 No 1–202 与原名（Species_original）。
  采纳名（Species_current_IF）单列给出——采纳名有 5 个重复，不能当主键。
* 删除依据 = 表2_裁决账本.xlsx「排除名单_29」的问题类别与排除原因。
* BUSCO 数值不做任何重算：按株一一对应搬运，只做删株、重编号、命名。
  唯一例外 = Umbelopsis_ramanniana 蛋白 complete_percent 一格（源表漏更新，项目已有裁定），见 OLD_FIX。

账目
----
  旧表 213 条  =  202 株现队列  +  11 株表2 已排除（非队列基因组）
  213 − 11 = 202

输入（全部只读）
----------------
1. 投稿附表/表3_研究总表.xlsx      → 权威 No / 原名 / 采纳名 / Subphylum / Assembly_accession
2. 投稿附表/表2_裁决账本.xlsx      → 排除名单_29 的问题类别与排除原因
3. protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv
   （213 长表最新副本，2026-08-17 22:53）→ BUSCO 数值；仅 1 处按项目既有裁定修正
   （Umbelopsis_ramanniana 蛋白 C 73.5 → 95.1，见 OLD_FIX 注释）

输出 → protein_BUSCO/comparison_busco_genome_longest_isoform_protein_202/
  longest_isoform_busco_genome_protein_mucoromycota_long_202.csv / .tsv   主表（202×8=1616 行）
  longest_isoform_busco_genome_protein_mucoromycota_wide_202.csv          宽表（202 行）
  longest_isoform_busco_genome_protein_mucoromycota_202.xlsx              长表+宽表+审计明细+说明
  审计变更明细_213_to_202.tsv                                             机器可读明细
  审计变更报告_213_to_202.md                                              人读报告
"""
import csv, re, os, collections
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ROOT = "/Users/caoyuan/Desktop/Genome_annotation_project_of_Mucoromycota"
T3 = f"{ROOT}/投稿附表/表3_研究总表.xlsx"
T2 = f"{ROOT}/投稿附表/表2_裁决账本.xlsx"
OLD = f"{ROOT}/protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv"
BA = os.path.expanduser("~/Desktop/omark_filtered_v6/omark_filtered_v6_qc/busco/busco_before_after.tsv")
OUTDIR = f"{ROOT}/protein_BUSCO/comparison_busco_genome_longest_isoform_protein_202"
os.makedirs(OUTDIR, exist_ok=True)
LONG_CSV = f"{OUTDIR}/longest_isoform_busco_genome_protein_mucoromycota_long_202.csv"
LONG_TSV = f"{OUTDIR}/longest_isoform_busco_genome_protein_mucoromycota_long_202.tsv"
WIDE_CSV = f"{OUTDIR}/longest_isoform_busco_genome_protein_mucoromycota_wide_202.csv"
DET_TSV = f"{OUTDIR}/审计变更明细_213_to_202.tsv"
REPORT = f"{OUTDIR}/审计变更报告_213_to_202.md"
XLSX = f"{OUTDIR}/longest_isoform_busco_genome_protein_mucoromycota_202.xlsx"


def collapse(name):
    """归一化：大小写/空格/下划线/点/连字符全部忽略。"""
    return re.sub(r'[\s_.\-]+', '', str(name)).lower()


# ------------------------------------------------------------------ 1. 权威骨干：表3
wb = openpyxl.load_workbook(T3, read_only=True)
ws = wb["研究菌株_202"]
backbone = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    backbone.append(dict(no=int(r[0]), orig=str(r[1]), strain=r[2],
                         accepted=str(r[3]), sub=r[4], acc=r[5], db=r[6]))
wb.close()
assert len(backbone) == 202 and [b['no'] for b in backbone] == list(range(1, 203))
print(f"[骨干] 表3 研究菌株_202：{len(backbone)} 株，No 连续 1–202")

# ------------------------------------------------------------------ 2. 排除账本：表2
wb = openpyxl.load_workbook(T2, read_only=True)
ws = wb["排除名单_29"]
excl_rows = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    excl_rows.append(dict(t1_no=r[0], orig=str(r[1]), accepted=r[2], acc=r[3], db=r[4],
                          biosample=r[5], C=r[6], D=r[7], category=r[8], reason=r[9],
                          keep=r[10], basis=r[11], note=r[12] if len(r) > 12 else None))
wb.close()
excl_by_orig = {}
for e in excl_rows:
    excl_by_orig.setdefault(collapse(e['orig']), []).append(e)
    excl_by_orig.setdefault(collapse(e['accepted']), []).append(e)
print(f"[账本] 表2 排除名单_29：{len(excl_rows)} 条")

# ------------------------------------------------------------------ 3. 旧表
# 源文件选择：protein_BUSCO/ 下有 3 份 213 长表副本——
#   · protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv       2026-08-17 22:53 ← 采用（最新）
#   · protein_BUSCO/final_comparison/…（同上，内容相同）
#   · protein_BUSCO/old_tables/comparison_busco_genome_longest_isoform_protein_213/…  2026-08-16 14:23（较旧）
# 新旧两份只差 1 株：Umbelopsis_ramanniana 的 longest-isoform 蛋白集（该株重新注释过）：
#   旧 72.9/0.6/9.2/17.3（C=73.5）  →  新 94.7/0.4/2.8/2.0
# 但新长表的 complete_percent 单元格漏更新（仍为 73.5，与 S+D=95.1 矛盾）。
# 该值项目内已有定论：protein_BUSCO/final_comparison/…_part4.csv 记 95.2，
# 原蛋白文件审查/Fig1b_说明与数据来源_v2.md（2026-09-19，问题1）裁定按 S+D 自洽值取 95.1，
# 且 Fig1b 派生表 braker_protein_BUSCO_C = 95.1 已据此出图。此处沿用 95.1，并在报告中显式记录。
OLD_FIX = {('Umbelopsis_ramanniana', 'P', 'complete_percent'): ('73.5', '95.1')}

old_rows = list(csv.DictReader(open(OLD, encoding='utf-8-sig')))
for r in old_rows:
    fx = OLD_FIX.get((r['sample'], r['dataset_label'], 'complete_percent'))
    if fx and r['complete_percent'] == fx[0]:
        r['complete_percent'] = fx[1]
old_samples = []
for r in old_rows:
    if r['sample'] not in old_samples:
        old_samples.append(r['sample'])
print(f"[旧表] {os.path.basename(OLD)}：{len(old_rows)} 行 / {len(old_samples)} 株")

by_sample = collections.defaultdict(list)
for i, r in enumerate(old_rows):
    r['_i'] = i
    by_sample[r['sample']].append(r)
# 保持旧表原始行序（每株 8 行：G-S, G-D, G-F, G-M, P-S, P-D, P-F, P-M），
# 不能用 stack_order 排序——那样会把 Genome 与 Protein 交错。
for s in by_sample:
    by_sample[s].sort(key=lambda r: r['_i'])

# ------------------------------------------------------------------ 4. 旧名 → No
# 旧 213 表的样本名混用了「原名」与「采纳名」（Mortierellaceae 的重新归类名如
# Actinomortierella / Gamsiella / Podila minutissima / Linnemannia 用的是采纳名）。
# 判定用「先到先得」四轮，保证每个 No 只被一条旧记录认领：
#   ① 命中表3 原名                     → 现队列该 No
#   ② 命中表3 采纳名 且该 No 尚未被认领 → 现队列该 No
#   ③ 命中表2 排除名单原名              → 表2 已排除（非队列基因组）
#   ④ 其余                              → 报错
# 关键例子 Linnemannia camargensis：它同时是 No69 的采纳名与表2 #209（JGI 镜像）的原名。
# 旧表里 Haplosporangium gracile（No69 原名）已在第①轮认领 No69，因此这条只能落到第③轮
# = 表2 已排除的镜像装配。这正是 213 = 202 + 11 的账。
by_orig, by_acc = {}, {}
for b in backbone:
    by_orig.setdefault(collapse(b['orig']), []).append(b)
    by_acc.setdefault(collapse(b['accepted']), []).append(b)
excl_by_orig, excl_by_acc = {}, {}
for e in excl_rows:
    excl_by_orig.setdefault(collapse(e['orig']), []).append(e)
    excl_by_acc.setdefault(collapse(e['accepted']), []).append(e)

assert not (set(by_orig) & set(excl_by_orig)), \
    f"表3 原名与表2 排除原名冲突：{set(by_orig) & set(excl_by_orig)}"

assign, excluded, match_kind, claimed = {}, [], {}, {}
ambiguous_log = []


def _pick(cands):
    return cands[0] if len(cands) == 1 else None


for s in old_samples:                                  # ① 原名
    b = _pick(by_orig.get(collapse(s), []))
    if b and b['no'] not in claimed:
        assign[s] = b; match_kind[s] = 'original'; claimed[b['no']] = s
for s in old_samples:                                  # ② 采纳名（No 未被认领）
    if s in assign:
        continue
    b = _pick(by_acc.get(collapse(s), []))
    if b and b['no'] not in claimed:
        assign[s] = b; match_kind[s] = 'accepted'; claimed[b['no']] = s
for s in old_samples:                                  # ③ 表2 排除名单
    if s in assign:
        continue
    if excl_by_orig.get(collapse(s)) or excl_by_acc.get(collapse(s)):
        excluded.append(s); match_kind[s] = 'excluded'
    else:
        raise SystemExit(f"无法判定旧表样本：{s}")

# 第②轮里 No 已被认领、只能在第③轮兜住的，单独记一笔（Linnemannia camargensis 就是这种）
for s in old_samples:
    if match_kind.get(s) == 'excluded':
        b = _pick(by_acc.get(collapse(s), []))
        if b:
            ambiguous_log.append((s, b, claimed.get(b['no'])))

unclaimed = [b['no'] for b in backbone if b['no'] not in claimed]
assert not unclaimed, f"这些 202 株没有认领到旧表记录：{unclaimed}"
print(f"[匹配] 命中现队列 {len(assign)}（覆盖 No {len(claimed)}/202）；已排除 {len(excluded)}")
for s, b, by in ambiguous_log:
    print(f"       · 同名辨析：{s} 虽为 No{b['no']} 的采纳名，但 No{b['no']} 已由 {by} 认领，"
          f"故判为表2 已排除的同名记录")

hits = collections.defaultdict(list)
for s, b in assign.items():
    hits[b['no']].append(s)
dupes = {no: ss for no, ss in hits.items() if len(ss) > 1}
assert not dupes, f"同一 No 被多个旧名命中：{dupes}"

# ------------------------------------------------------------------ 5. 组装 202 株长表
new_long, audit = [], []
for b in backbone:
    no = b['no']
    cands = hits.get(no, [])
    assert cands, f"No{no} 在旧表中没有对应样本"
    keep = cands[0]
    rows = by_sample[keep]
    assert len(rows) == 8, (keep, len(rows))
    assert [r['metric'] for r in rows] == ['S', 'D', 'F', 'M'] * 2, (keep, [r['metric'] for r in rows])
    assert [r['dataset_label'] for r in rows] == ['G'] * 4 + ['P'] * 4, keep
    renamed = collapse(b['orig']) != collapse(b['accepted'])
    name_status = 'renamed_in_current_taxonomy' if renamed else 'unchanged'
    for r in rows:
        new_long.append({
            'sample': b['orig'],                 # 主名 = 原名（同表3/表5/表6/表7）
            'sample_order': no,                  # 表3 权威顺序 1–202
            'lineage': r['lineage'],
            'dataset': r['dataset'],
            'dataset_label': r['dataset_label'],
            'metric': r['metric'],
            'category': r['category'],
            'stack_order': r['stack_order'],
            'percent': r['percent'],
            'complete_percent': r['complete_percent'],
            'n_buscos': r['n_buscos'],
            'no': no,
            'species_original': b['orig'],
            'species_current_IF': b['accepted'],
            'subphylum': b['sub'],
            'assembly_accession': b['acc'],
            'name_status': name_status,
            'source_sample_in_213': keep,
        })
    audit.append(dict(no=no, species_original=b['orig'], species_current_IF=b['accepted'],
                      subphylum=b['sub'], assembly_accession=b['acc'],
                      status='kept', name_status=name_status,
                      source_sample_in_213=keep, exclusion_category='', exclusion_reason='',
                      kept_representative='', note=''))

# 被排除的 11 条
del_detail = []
for s in sorted(excluded):
    k = collapse(s)
    e = (excl_by_orig.get(k) or excl_by_acc.get(k) or [None])[0]
    cat = e['category'] if e else ''
    reason = e['reason'] if e else '不在表3 研究菌株_202 中'
    keeprep = e['keep'] if e else ''
    t1 = e['t1_no'] if e else ''
    # 排除账本里还有别的记录也以该名字作为采纳名（如 No69），提示这条是镜像/同种的那一份
    twin = [b for b in backbone if collapse(b['accepted']) == k]
    note = ''
    if twin:
        note = (f"该条旧表用的是采纳名，实为表2 编号 {t1} 的 {e['orig'] if e else ''}"
                f"（{cat}，{keeprep and '保留代表 ' + str(keeprep) or '不保留'}）；"
                f"与现队列 No{twin[0]['no']}（{twin[0]['orig']}）为同一物种的镜像装配，"
                f"新表只保留 No{twin[0]['no']} 那一份（NCBI 侧）")
    del_detail.append(dict(old_sample=s, t1_no=t1, orig=e['orig'] if e else s.replace('_', ' '),
                           accepted=e['accepted'] if e else '', acc=e['acc'] if e else '',
                           db=e['db'] if e else '', C=e['C'] if e else '', D=e['D'] if e else '',
                           category=cat, reason=reason, keep=keeprep, note=note))
    audit.append(dict(no='', species_original=e['orig'] if e else s.replace('_', ' '),
                      species_current_IF=e['accepted'] if e else '',
                      subphylum='', assembly_accession=e['acc'] if e else '',
                      status='DELETE', name_status='excluded',
                      source_sample_in_213=s, exclusion_category=cat,
                      exclusion_reason=reason, kept_representative=keeprep, note=note))

assert len(new_long) == 202 * 8 and len(del_detail) == 11, (len(new_long), len(del_detail))
print(f"[账目] 213 = {len(hits)} 现队列 + {len(del_detail)} 已排除 → 新表 {len(hits)} 株")

# ------------------------------------------------------------------ 6. 校验
assert len(set(r['sample'] for r in new_long)) == 202
assert set(collections.Counter(r['sample'] for r in new_long).values()) == {8}
assert set(r['lineage'] for r in new_long) == {'mucoromycota'}
assert set(r['n_buscos'] for r in new_long) == {'1560'}
grp = collections.defaultdict(float)
for r in new_long:
    grp[(r['sample'], r['dataset'])] += float(r['percent'])
bad = {k: v for k, v in grp.items() if abs(v - 100) > 0.15}
assert not bad, bad
# C 必须与 S+D 自洽（允许 0.15 的四舍五入差）
csd = []
for (s, d) in {(r['sample'], r['dataset']) for r in new_long}:
    q = {(r['metric']): r for r in new_long if r['sample'] == s and r['dataset'] == d}
    c = float(q['S']['complete_percent'])
    if abs(c - (float(q['S']['percent']) + float(q['D']['percent']))) > 0.15:
        csd.append((s, d, c, float(q['S']['percent']) + float(q['D']['percent'])))
assert not csd, csd
print("[校验] 每株 C = S + D（≤0.15 舍入差）✓")
for r in new_long:
    src = [x for x in by_sample[r['source_sample_in_213']]
           if x['dataset'] == r['dataset'] and x['metric'] == r['metric']]
    assert len(src) == 1
    assert (src[0]['percent'], src[0]['complete_percent'], src[0]['n_buscos']) == \
           (r['percent'], r['complete_percent'], r['n_buscos']), (r, src[0])
# 位置级逐格比对：每株的 8 行必须与旧表那一株的 8 行按同序、同值拷贝
KEEPCOLS = ['sample', 'lineage', 'dataset', 'dataset_label', 'metric', 'category',
            'stack_order', 'percent', 'complete_percent', 'n_buscos']
for b in backbone:
    new8 = [r for r in new_long if r['no'] == b['no']]
    old8 = by_sample[new8[0]['source_sample_in_213']]
    assert [tuple(r[c] for c in KEEPCOLS[1:]) for r in new8] == \
           [tuple(r[c] for c in KEEPCOLS[1:]) for r in old8], b['no']
print("[校验] 每株 8 行与旧表按同序逐格相同（含 dataset/metric/stack_order/percent/C/n）✓")

# 独立交叉校验：基因组 C 必须等于表3 自带的 BUSCO C_muco 列（该列独立记录于表3，非本表来源）
wb = openpyxl.load_workbook(T3, read_only=True)
ws = wb["研究菌株_202"]
c_muco = {}
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    c_muco[int(r[0])] = (str(r[1]), r[8])
wb.close()
bad_c = []
for b in backbone:
    orig, c = c_muco[b['no']]
    row = [r for r in new_long if r['no'] == b['no'] and r['dataset_label'] == 'G'][0]
    if orig != row['sample']:
        bad_c.append((b['no'], '名不符', orig, row['sample']))
    elif c is not None and abs(float(c) - float(row['complete_percent'])) > 0.051:
        bad_c.append((b['no'], 'C 不符', c, row['complete_percent']))
assert not bad_c, bad_c
print(f"[校验] 基因组 C 与表3 独立记录的 BUSCO C_muco 全部一致（{sum(1 for v in c_muco.values() if v[1] is not None)}/202 有值）✓")
print("[校验] 202 株 × 8 行 = 1616；每株 8 行；S+D+F+M=100；数值与旧表逐格一致 ✓")

# ------------------------------------------------------------------ 7. 写长表
COLS = ['sample', 'sample_order', 'lineage', 'dataset', 'dataset_label', 'metric',
        'category', 'stack_order', 'percent', 'complete_percent', 'n_buscos',
        'no', 'species_original', 'species_current_IF', 'subphylum',
        'assembly_accession', 'name_status', 'source_sample_in_213']
with open(LONG_CSV, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=COLS); w.writeheader(); w.writerows(new_long)
with open(LONG_TSV, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=COLS, delimiter='\t'); w.writeheader(); w.writerows(new_long)
print("[写出]", os.path.basename(LONG_CSV), "/", os.path.basename(LONG_TSV))

# ------------------------------------------------------------------ 8. 宽表
wide_rows = []
for b in backbone:
    rows = [r for r in new_long if r['no'] == b['no'] and r['dataset_label'] == 'G']
    prows = [r for r in new_long if r['no'] == b['no'] and r['dataset_label'] == 'P']
    assert len(rows) == 4 and len(prows) == 4
    gg = {r['metric']: r['percent'] for r in rows}
    pp = {r['metric']: r['percent'] for r in prows}
    wide_rows.append({
        'no': b['no'], 'sample': b['orig'], 'species_current_IF': b['accepted'],
        'subphylum': b['sub'], 'assembly_accession': b['acc'],
        'name_status': 'renamed_in_current_taxonomy'
                       if collapse(b['orig']) != collapse(b['accepted']) else 'unchanged',
        # C 直接沿用长表 complete_percent：BUSCO 的 C 由计数算出后各自四舍五入，
        # 少数株 C 与 S+D 差 0.1（如 Absidia healeyae 蛋白 C=94.4 vs S+D=94.5）。
        'C_genome': rows[0]['complete_percent'],
        'S_genome': gg['S'], 'D_genome': gg['D'], 'F_genome': gg['F'], 'M_genome': gg['M'],
        'C_protein': prows[0]['complete_percent'],
        'S_protein': pp['S'], 'D_protein': pp['D'], 'F_protein': pp['F'], 'M_protein': pp['M'],
        'delta_C_protein_minus_genome': str(round(float(prows[0]['complete_percent'])
                                                  - float(rows[0]['complete_percent']), 1)),
        'n_buscos': rows[0]['n_buscos'],
    })
with open(WIDE_CSV, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(wide_rows[0].keys())); w.writeheader(); w.writerows(wide_rows)
print("[写出]", os.path.basename(WIDE_CSV))
for wr in wide_rows:
    s = [r for r in new_long if r['no'] == wr['no'] and r['dataset_label'] == 'G'][0]
    assert wr['C_genome'] == s['complete_percent']
    s = [r for r in new_long if r['no'] == wr['no'] and r['dataset_label'] == 'P'][0]
    assert wr['C_protein'] == s['complete_percent']
print("[校验] 宽表 C = 长表 complete_percent ✓")

# ------------------------------------------------------------------ 9. 审计明细 TSV
ACOLS = ['no', 'species_original', 'species_current_IF', 'subphylum', 'assembly_accession',
         'status', 'name_status', 'source_sample_in_213',
         'exclusion_category', 'exclusion_reason', 'kept_representative', 'note']
audit.sort(key=lambda a: (1 if a['status'] == 'DELETE' else 0, a['no'] if a['no'] != '' else 0))
with open(DET_TSV, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=ACOLS, delimiter='\t'); w.writeheader(); w.writerows(audit)
print("[写出]", os.path.basename(DET_TSV))

# ------------------------------------------------------------------ 10. 报告 md
renamed = [b for b in backbone if collapse(b['orig']) != collapse(b['accepted'])]
cat_count = collections.Counter(d['category'] for d in del_detail)
lines = []
A = lines.append

A("# 213 → 202：BUSCO 基因组–蛋白对比表（mucoromycota_odb12）审计变更报告")
A("")
A("本报告说明 `longest_isoform_busco_genome_protein_mucoromycota_long.csv`（旧表，213 株）")
A("如何按当前审计重建为 202 株的新表。**旧表未做任何改动**（只读拷贝）。")
A("")
A("## 0. 数据来源与 2 处历史遗留")
A("")
A("### 0.1 用了哪一份 213 长表")
A("")
A("`protein_BUSCO/` 下有 3 份内容近乎相同的 213 长表副本，本表用的是**最新的一份**：")
A("")
A("| 副本 | 时间 | 是否采用 |")
A("|---|---|---|")
A("| `protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv` | 2026-08-17 22:53 | ✅ 采用 |")
A("| `protein_BUSCO/final_comparison/` 同名文件 | 2026-08-17 22:53 | 内容相同 |")
A("| `protein_BUSCO/old_tables/comparison_busco_genome_longest_isoform_protein_213/` 同名文件 | 2026-08-16 14:23 | ❌ 较旧，未采用 |")
A("")
A("新旧两份**只差 1 株**——Umbelopsis ramanniana 的 longest-isoform 蛋白集（该株做过重新注释）：")
A("")
A("```")
A("旧副本  蛋白 S/D/F/M = 72.9 / 0.6 / 9.2 / 17.3   C = 73.5")
A("新副本  蛋白 S/D/F/M = 94.7 / 0.4 / 2.8 /  2.0   C = 73.5  ← C 单元格漏改")
A("```")
A("")
A("### 0.2 已按项目既有裁定修正 1 个单元格")
A("")
A("新副本里 Umbelopsis ramanniana 蛋白的 `complete_percent` 仍是旧的 73.5，与同行 S+D=95.1 矛盾。")
A("该问题项目内已处理过：")
A("")
A("- `protein_BUSCO/final_comparison/…_part4.csv`（绘图用 4 分表）记该株蛋白 **C = 95.2**；")
A("- `原蛋白文件审查/Fig1b_说明与数据来源_v2.md`（2026-09-19「问题 1（已处理）」）裁定按 S+D 自洽值取 **95.1**，")
A("  且 Fig1b 派生表 `braker_protein_BUSCO_C = 95.1` 已据此出图。")
A("")
A("本表沿用项目已出图的口径 **95.1**（与 parts 的 95.2 差 0.1，属 BUSCO C 与 S+D 各自四舍五入）。")
A("除这一格外，全表数值与源表逐格相同。")
A("")
A("## 1. 一句话账目")
A("")
A("```")
A(f"旧表          {len(old_samples)} 株  × 2 数据集(Genome / Longest-isoform protein) × 4 类别 = {len(old_rows)} 行")
A(f"删（排除）    -{len(del_detail)} 株（表2 裁决账本 排除名单_29 之内的非队列基因组）")
A(f"保留          {len(hits)} 株  × 8 = {len(new_long)} 行   ← 即当前 202 株标准")
A(f"改名          {len(renamed)} 株（原名 → 采纳名，Index Fungorum 现名）")
A("```")
A("")
A("## 2. 删除的 11 株（旧表有、202 标准没有）")
A("")
A("| 旧表样本名 | 表1编号 | 采纳名 | 问题类别 | 排除原因 | 保留代表 |")
A("|---|---|---|---|---|---|")
for d in del_detail:
    A(f"| `{d['old_sample']}` | {d['t1_no']} | {d['accepted'] or '—'} | {d['category']} | "
      f"{str(d['reason']).replace('|', '/')} | {d['keep'] or '—'} |")
A("")
A("类别构成：" + "、".join(f"{k} {v}" for k, v in sorted(cat_count.items())) + "。")
A("")
for d in del_detail:
    if d['note']:
        A(f"> **注意 `{d['old_sample']}`**：{d['note']}。")
        A("")
A("## 3. 改名的 15 株（原名 → 采纳名）")
A("")
A("新表的 `sample` 列沿用**原名**（与表3/表5/表6/表7 主名口径一致、且原名唯一），")
A("采纳名放在 `species_current_IF` 列；`name_status=renamed_in_current_taxonomy` 标记。")
A("")
A("| No | 原名（新表 sample） | 采纳名（species_current_IF） |")
A("|---|---|---|")
for b in renamed:
    A(f"| {b['no']} | {b['orig']} | {b['accepted']} |")
A("")
A("其中 5 个采纳名在 202 株内部**重名**（各自对应两株不同的菌株，因此采纳名不能当主键）：")
A("")
dup_acc = {k: [b['no'] for b in backbone if b['accepted'] == k]
           for k, v in collections.Counter(b['accepted'] for b in backbone).items() if v > 1}
for k, nos in sorted(dup_acc.items()):
    det = "；".join(f"No{n} {[b['orig'] for b in backbone if b['no'] == n][0]}" for n in nos)
    A(f"- **{k}** ← {det}")
A("")
A("## 4. 新表怎么用")
A("")
A(f"- 主表：`{os.path.basename(LONG_CSV)}` / `.tsv`，1616 行 = 202 株 × 8；")
A("  列 = 旧表原有列（`sample, sample_order, lineage, dataset, dataset_label, metric, category, "
  "stack_order, percent, complete_percent, n_buscos`）+ 审计列（`no, species_original, "
  "species_current_IF, subphylum, assembly_accession, name_status, source_sample_in_213`）。")
A(f"- 宽表：`{os.path.basename(WIDE_CSV)}`，202 行，含 `delta_C_protein_minus_genome`。")
A(f"- Excel：`{os.path.basename(XLSX)}`（长表 / 宽表 / 审计明细 / 说明 四个 sheet）。")
A(f"- 机器可读明细：`{os.path.basename(DET_TSV)}`。")
A("- `sample_order` 已按表3 权威顺序重编号为 **1–202**（旧表为 1–213）；")
A("  绘图脚本 `plot_longest_isoform_busco_original_order_v8_underscore_fix.r` 按行序取样本顺序，会自动跟随。")
A("- BUSCO 数值**未做任何重算**，逐格等于旧表；仅 202 株的子集。")
A("- 独立交叉校验：新表 202 株的基因组 C 与表3 自带的 `BUSCO C_muco` 列（独立记录）**全部一致**，")
A("  说明每株都取到了正确基因组的数值。")
A("")

# OMArk 附注
if os.path.exists(BA):
    ba = list(csv.DictReader(open(BA, encoding='utf-8'), delimiter='\t'))
    A("## 5. 待你确认：10 株 OMArk 过滤后的蛋白 BUSCO 未并入")
    A("")
    A("2026-10-05 对以下 10 株做了 OMArk 污染过滤后的蛋白集重跑（仅蛋白，基因组未重跑），")
    A("数值与旧表不同。本表为保持全表口径统一，**沿用旧表数值**；若要改用过滤后数值，请说明，另出一版。")
    A("")
    A("| 样本 | 旧表/过滤前 C | 过滤后 C | 过滤前 S/D | 过滤后 S/D |")
    A("|---|---|---|---|---|")
    for r in ba:
        A(f"| {r['sample']} | {r['before_C']} | {r['after_C']} | {r['before_S']}/{r['before_D']} | "
          f"{r['after_S']}/{r['after_D']} |")
    A("")
    A(f"来源：`{BA}`")
    A("")

A("## 6. 复现")
A("")
A("```bash")
A(f"python3 {os.path.relpath(os.path.abspath(__file__), ROOT)}")
A("```")
A("")
A("输入（只读）：表3_研究总表.xlsx `研究菌株_202`、表2_裁决账本.xlsx `排除名单_29`、旧 213 长表。")
open(REPORT, 'w', encoding='utf-8').write("\n".join(lines) + "\n")
print("[写出]", os.path.basename(REPORT))

# ------------------------------------------------------------------ 11. xlsx
wb = openpyxl.Workbook()
hdr_font = Font(bold=True, color="FFFFFF", size=10)
hdr_fill = PatternFill("solid", fgColor="1F4E79")
thin = Side(style="thin", color="BFBFBF")
border = Border(left=thin, right=thin, top=thin, bottom=thin)
warn_fill = PatternFill("solid", fgColor="FFF2CC")
ren_fill = PatternFill("solid", fgColor="E2EFDA")

ws1 = wb.active; ws1.title = "1_长表_202株"
ws1.append(COLS)
for r in new_long:
    ws1.append([r[c] for c in COLS])
ws1.freeze_panes = "A2"
ws1.auto_filter.ref = f"A1:{get_column_letter(len(COLS))}{ws1.max_row}"
for i, wd in enumerate([24, 8, 12, 9, 22, 7, 22, 7, 8, 10, 9, 5, 24, 24, 18, 18, 26, 24], 1):
    ws1.column_dimensions[get_column_letter(i)].width = wd
    c = ws1.cell(row=1, column=i); c.font = hdr_font; c.fill = hdr_fill
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True); c.border = border

ws2 = wb.create_sheet("2_宽表_202株")
wcols = list(wide_rows[0].keys())
ws2.append(wcols)
for r in wide_rows:
    ws2.append([r[c] for c in wcols])
ws2.freeze_panes = "A2"
for i in range(1, len(wcols) + 1):
    ws2.column_dimensions[get_column_letter(i)].width = 22
    c = ws2.cell(row=1, column=i); c.font = hdr_font; c.fill = hdr_fill
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True); c.border = border

ws3 = wb.create_sheet("3_审计变更明细")
ws3.append(ACOLS)
for a in audit:
    ws3.append([a[c] for c in ACOLS])
ws3.freeze_panes = "A2"
for i in range(1, len(ACOLS) + 1):
    ws3.column_dimensions[get_column_letter(i)].width = 26 if i in (2, 3, 10, 12) else 18
    c = ws3.cell(row=1, column=i); c.font = hdr_font; c.fill = hdr_fill
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True); c.border = border
for row in ws3.iter_rows(min_row=2, max_row=ws3.max_row, max_col=len(ACOLS)):
    if row[5].value == 'DELETE':
        for c in row:
            c.fill = warn_fill
    elif row[6].value == 'renamed_in_current_taxonomy':
        for c in row:
            c.fill = ren_fill

ws4 = wb.create_sheet("4_说明")
notes = [
    ("202 株 BUSCO 基因组–蛋白对比表（mucoromycota_odb12）— 说明", True),
    ("", False),
    ("• 数据来源：protein_BUSCO/longest_isoform_busco_genome_protein_mucoromycota_long.csv"
     "（213 株长表最新副本，2026-08-17 22:53；只读，未改动）。注意 protein_BUSCO/old_tables/… 下那份是"
     "2026-08-16 的较旧副本（Umbelopsis ramanniana 蛋白仍是重新注释前的值），不要用。", False),
    ("• 本表 = 旧表按当前审计裁剪并重编号后的 202 株版本；除下条注明的 1 格外，BUSCO 数值一律照搬，未做任何重算。", False),
    ("• 已修正 1 格：Umbelopsis ramanniana 蛋白 complete_percent 73.5 → 95.1。源表该格漏更新（与同行 S+D=95.1 矛盾）；"
     "项目既有裁定见 原蛋白文件审查/Fig1b_说明与数据来源_v2.md「问题1（已处理）」，"
     "绘图用 part4 分表记 95.2（差 0.1 属 BUSCO C 与 S+D 各自四舍五入）。", False),
    ("• 213 = 202（现队列）+ 11（表2 裁决账本已排除的非队列基因组）。11 株明细见「3_审计变更明细」"
     "（黄色行）。其中 Linnemannia_camargensis 是 No69 的 JGI 镜像装配（MIRROR_COPY），"
     "在旧表里以采纳名出现，故看似「同一株跑了两遍」；新表只保留 No69（Haplosporangium gracile）。", False),
    ("• 主键：No（1–202，表3 权威编号）+ Species_original（原名，202 个互不重复）。采纳名有 5 个重复"
     "（Cunninghamella elegans / Linnemannia hyalina / Mucor circinelloides / Rhizopus arrhizus / "
     "Rhizopus microsporus），不能作主键，故单列为 Species_current_IF。", False),
    ("• sample 列 = Species_original（原名），与表3/表5/表6/表7 主名口径一致；source_sample_in_213 "
     "记录该行在旧表中的样本名。", False),
    ("• 改名 15 株：name_status = renamed_in_current_taxonomy（绿色行），采纳名见 species_current_IF。", False),
    ("• sample_order 已按表3 权威顺序重编号为 1–202（旧表为 1–213）。", False),
    ("• lineage 全部为 mucoromycota（odb12，n=1560）；dataset=Genome / Longest-isoform protein，"
     "dataset_label=G / P；metric：S 单拷贝完整、D 重复完整、F 片段化、M 缺失；complete_percent = C = S + D。", False),
    ("• 宽表的 C 直接取长表 complete_percent（个别株 C 与 S+D 差 0.1，属 BUSCO 各自四舍五入）。", False),
    ("", False),
    ("⚠ 需要留意（本表未并入）：", True),
    ("• 2026-10-05 对 10 株 OMArk 污染过滤后的蛋白集（Desktop/omark_filtered_v6/）重跑了蛋白 BUSCO，"
     "数值与旧表不同（见 omark_filtered_v6_qc/busco/busco_before_after.tsv）。本表沿用旧表数值以保证"
     "全表口径统一；若要改用过滤后数值，请说明，另出一版。", False),
]
for i, (text, bold) in enumerate(notes, 1):
    c = ws4.cell(row=i, column=1, value=text)
    c.font = Font(bold=bold, size=10)
    c.alignment = Alignment(wrap_text=True, vertical="top")
ws4.column_dimensions['A'].width = 150

wb.save(XLSX)
print("[写出]", os.path.basename(XLSX))

# ------------------------------------------------------------------ 12. 控制台摘要
print("\n===== 删除 %d 株（表2 排除名单内的非队列基因组）=====" % len(del_detail))
for d in del_detail:
    print(f"   - {d['old_sample']:<34} [{d['category']}] {str(d['reason'])[:56]}")
print("\n===== 改名 %d 株 =====" % len(renamed))
for b in renamed:
    print(f"   - No{b['no']:>3}  {b['orig']}  ->  {b['accepted']}")
print("\n完成。输出目录:", OUTDIR)
