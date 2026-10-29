#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""10_OMARK/omark_final_per_strain.tsv（212 株）→ 202 株审计版

原始表格只读，不动一个字节；本脚本读它、按株拷贝、重命名后另存新文件。

212 是怎么来的
--------------
  212 行 = 202（现队列） + 9（表2 已排除、且当时进了 OMArk 的株） + 1（No69 的镜像重复条目）
  202 行 = 192（用 run=original 行） + 10（用 run=filtered 行，即 OMArk 污染过滤后重跑的那 10 株）
  212 − 9 − 1 = 202

处理三件事
----------
1. **选行**：同一株只用一行。10 株过滤过的用 run=filtered 那行（这是当前结果），
   其余 192 株用 run=original 行。全表 is_final 均为 True。
2. **删行**：9 株表2 已排除 + 1 条 No69 镜像重复（Linnemannia camargensis，
   采纳名与 No69 同名，实为表2 #209 的 JGI 镜像装配）。
3. **命名**：第 1 列 sample 保持原始表写法（= 表3 原名，202 个唯一）；
   第 2 列新增 `现用名` = 采纳名，但下列 6 株因采纳名与表内另一株重名而**保留原名**：
       No41 Cunninghamella bertholletiae、No94 Mortierella hygrophila、
       No111 Mucor ambiguus、No127 Mucor velutinosus、
       No151 Rhizopus azygosporus、No152 Rhizopus delemar。
   其余 9 株改名生效。这样 202 个现用名互不重复。

输出 → 10_OMARK/
  omark_final_per_strain_202.tsv / .csv / .xlsx
  omark_202_处理说明.md
"""
import csv, re, os, collections
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ROOT = "/Users/caoyuan/Desktop/Genome_annotation_project_of_Mucoromycota"
SRC = f"{ROOT}/10_OMARK/omark_final_per_strain.tsv"
T3 = f"{ROOT}/投稿附表/表3_研究总表.xlsx"
T2 = f"{ROOT}/投稿附表/表2_裁决账本.xlsx"
OUTDIR = f"{ROOT}/10_OMARK"
OUT_TSV = f"{OUTDIR}/omark_final_per_strain_202.tsv"
OUT_CSV = f"{OUTDIR}/omark_final_per_strain_202.csv"
OUT_XLSX = f"{OUTDIR}/omark_final_per_strain_202.xlsx"
REPORT = f"{OUTDIR}/omark_202_处理说明.md"

REVERT_TO_ORIGINAL = {41, 94, 111, 127, 151, 152}


def collapse(n):
    return re.sub(r'[\s_.\-]+', '', str(n)).lower()


# ================================================================== 1. 读原始表（只读）
rows = list(csv.DictReader(open(SRC, encoding='utf-8'), delimiter='\t'))
ORIG_COLS = list(rows[0].keys())
n_uniq_src = len(set(r['sample'] for r in rows))
print(f"[原始表] {os.path.basename(SRC)}：{len(rows)} 行 × {len(ORIG_COLS)} 列（只读）")
run_cnt = collections.Counter(r['run'] for r in rows)
print(f"[原始表] run 分布：{dict(run_cnt)}")

# ================================================================== 2. 表3 骨干 + 表2 排除账本
wb = openpyxl.load_workbook(T3, read_only=True)
ws = wb["研究菌株_202"]
backbone = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    backbone.append(dict(no=int(r[0]), orig=str(r[1]), accepted=str(r[3]),
                         sub=r[4], acc=r[5]))
wb.close()
assert len(backbone) == 202 and [b['no'] for b in backbone] == list(range(1, 203))
no2b = {b['no']: b for b in backbone}

wb = openpyxl.load_workbook(T2, read_only=True)
ws = wb["排除名单_29"]
excl_rows = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    excl_rows.append(dict(t1_no=r[0], orig=str(r[1]), accepted=r[2], acc=r[3],
                          category=r[8], reason=r[9], keep=r[10]))
wb.close()
excl_by_orig, excl_by_acc = {}, {}
for e in excl_rows:
    excl_by_orig.setdefault(collapse(e['orig']), []).append(e)
    excl_by_acc.setdefault(collapse(e['accepted']), []).append(e)
print(f"[骨干] 表3 202 株；[账本] 表2 排除名单_29 {len(excl_rows)} 条")

by_orig, by_acc = {}, {}
for b in backbone:
    by_orig.setdefault(collapse(b['orig']), []).append(b)
    by_acc.setdefault(collapse(b['accepted']), []).append(b)

# 现用名：采纳名，6 株重名的保留原名
cur_of = {b['no']: (b['orig'] if b['no'] in REVERT_TO_ORIGINAL else b['accepted'])
          for b in backbone}
assert len(set(cur_of.values())) == 202
assert len(set(b['orig'] for b in backbone)) == 202
n_diff = sum(1 for b in backbone if cur_of[b['no']] != b['orig'])

# ================================================================== 3. 样本身份判定
#   ① 命中表3 原名             → 该 No
#   ② 命中表3 采纳名（No 未被认领）→ 该 No
#   ③ 命中表2 排除名单          → 删除
# ② 必须带「No 未被认领」条件：Linnemannia camargensis 既是 No69 的采纳名，
#    又是表2 #209（JGI 镜像）的原名；旧表里 Haplosporangium gracile（No69 原名）
#    已在第①轮认领 No69，所以这一条只能落到第③轮 = 删除。
claimed, dropped = {}, []
for r in rows:                                            # ① 原名
    k = collapse(r['sample'])
    c = by_orig.get(k, [])
    if len(c) == 1 and c[0]['no'] not in claimed:
        claimed[c[0]['no']] = r
for r in rows:                                            # ② 采纳名
    k = collapse(r['sample'])
    if any(x['sample'] == r['sample'] for x in claimed.values()):
        continue
    c = by_acc.get(k, [])
    if len(c) == 1 and c[0]['no'] not in claimed:
        claimed[c[0]['no']] = r
for r in rows:                                            # ③ 表2 排除
    if any(x['sample'] == r['sample'] for x in claimed.values()):
        continue
    k = collapse(r['sample'])
    e = (excl_by_orig.get(k) or excl_by_acc.get(k) or [None])[0]
    assert e, f"无法判定：{r['sample']}"
    dropped.append((r, e))
missing = [b['no'] for b in backbone if b['no'] not in claimed]
assert not missing, f"这些 202 株没有拿到 OMArk 结果：{missing}"
assert len(dropped) == 10, len(dropped)
print(f"[判定] 202 株全部拿到结果；删除 {len(dropped)} 行 → 212 − {len(dropped)} = {len(rows)-len(dropped)}")
print(f"[选行] run 分布（保留行）："
      f"{dict(collections.Counter(claimed[b['no']]['run'] for b in backbone))}")

# ================================================================== 4. 组装 202 行
OUT_COLS = ['sample', '现用名'] + [c for c in ORIG_COLS if c != 'sample']
# 第 1 列 sample = 原表原样照搬（下划线风格，保证可回溯到 OMArk 的 sum_file/run）
# 第 2 列 现用名 = 采纳名，6 株重名的保留原名；用与原表一致的「空格→下划线」写法
# 现用名的写法：
#   · 保持原名的株（187 未改 + 6 株重名保留）→ 直接沿用原表写法，避免出现「假改名」
#     （原表有 6 处下划线写法与表3 不同，如 Absidia_sp.YW-2025a ↔ Absidia sp. YW-2025a，
#      这只是分隔符差异，不是改名）
#   · 真正换名的株 → 用采纳名，并按原表风格把空格写成下划线
out_rows, out_by_no, disp_style = [], {}, []
for b in backbone:
    src = claimed[b['no']]
    if b['no'] in REVERT_TO_ORIGINAL or collapse(src['sample']) == collapse(b['accepted']):
        cur_disp = src['sample']
    else:
        cur_disp = b['accepted'].replace(' ', '_')
        assert cur_disp != src['sample'], b['no']
    assert collapse(cur_disp) == collapse(cur_of[b['no']]), b['no']
    if collapse(src['sample']) != collapse(b['orig']):
        disp_style.append((b['no'], b['orig'], src['sample']))
    o = {'sample': src['sample'], '现用名': cur_disp}
    for c in ORIG_COLS:
        if c != 'sample':
            o[c] = src[c]
    out_rows.append(o)
    out_by_no[b['no']] = o
n_real_rename = sum(1 for r in out_rows if r['sample'] != r['现用名'])
print(f"[写法] 原表与表3 写法不同的株（仅分隔符差异）：{len(disp_style)}")
for no, orig, srcname in disp_style:
    print(f"        No{no}: 表3 '{orig}'  原表 '{srcname}'")
print(f"[命名] 现用名 与 sample 实际不同的：{n_real_rename} 株")
assert len(out_rows) == 202
assert len(set(r['sample'] for r in out_rows)) == 202
assert len(set(r['现用名'] for r in out_rows)) == 202

# ================================================================== 5. 写表
with open(OUT_TSV, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter='\t'); w.writeheader(); w.writerows(out_rows)
with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=OUT_COLS); w.writeheader(); w.writerows(out_rows)
print(f"[写出] {os.path.basename(OUT_TSV)} / .csv（{len(OUT_COLS)} 列 × {len(out_rows)} 行）")

# xlsx
hdr_font = Font(bold=True, color="FFFFFF", size=10)
hdr_fill = PatternFill("solid", fgColor="1F4E79")
thin = Side(style="thin", color="BFBFBF")
border = Border(left=thin, right=thin, top=thin, bottom=thin)
hl = PatternFill("solid", fgColor="FFF2CC")
wb = openpyxl.Workbook()
ws1 = wb.active; ws1.title = "1_OMArk_202株"
ws1.append(OUT_COLS)
for r in out_rows:
    ws1.append([r[c] for c in OUT_COLS])
ws1.freeze_panes = "C2"
ws1.auto_filter.ref = f"A1:{get_column_letter(len(OUT_COLS))}{ws1.max_row}"
for i in range(1, len(OUT_COLS) + 1):
    ws1.column_dimensions[get_column_letter(i)].width = 26 if i <= 3 else 14
    c = ws1.cell(row=1, column=i); c.font = hdr_font; c.fill = hdr_fill
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    c.border = border
for row in ws1.iter_rows(min_row=2, max_row=ws1.max_row, max_col=len(OUT_COLS)):
    if row[2].value == 'filtered':           # run 列
        for c in row:
            c.fill = hl

ws2 = wb.create_sheet("2_处理说明")
notes = [
    ("OMArk 每株结果表 · 202 株审计版 — 说明", True),
    ("", False),
    (f"• 来源：10_OMARK/omark_final_per_strain.tsv（{len(rows)} 行 × {len(ORIG_COLS)} 列，只读，未改动）。", False),
    ("• 212 = 202（现队列）+ 9（表2 已排除且当时进了 OMArk 的株）+ 1（No69 的镜像重复条目）；"
     "212 − 10 = 202。", False),
    ("• 选行：10 株用 run=filtered 行（OMArk 污染过滤后重跑，黄色行）；其余 192 株用 run=original 行。", False),
    ("• 命名：第 1 列 sample = 原名（同表3 原名，202 个唯一）；第 2 列 现用名 = 采纳名，"
     "但 6 株因采纳名与表内另一株重名而保留原名（No41/94/111/127/151/152）。"
     "两列只在 9 株上不同。", False),
    ("• 除新增的「现用名」列外，其余 50 列一律原样照搬，未做任何换算。", False),
    ("• 删除的 10 行见 omark_202_处理说明.md。", False),
]
for i, (t, bold) in enumerate(notes, 1):
    c = ws2.cell(row=i, column=1, value=t)
    c.font = Font(bold=bold, size=10)
    c.alignment = Alignment(wrap_text=True, vertical="top")
ws2.column_dimensions['A'].width = 150
wb.save(OUT_XLSX)
print(f"[写出] {os.path.basename(OUT_XLSX)}")

# ================================================================== 6. 报告
L = []
A = L.append
A("# 10_OMARK/omark_final_per_strain.tsv → 202 株审计版 · 处理说明")
A("")
A(f"原始表：`10_OMARK/omark_final_per_strain.tsv`（{len(rows)} 行 × {len(ORIG_COLS)} 列）——**只读，未改动**。")
A(f"新表：`10_OMARK/omark_final_per_strain_202.tsv` / `.csv` / `.xlsx`（**{len(out_rows)} 行 × {len(OUT_COLS)} 列**）。")
A("")
A("## 0. 一句话账目")
A("")
A("```")
A(f"原始表            {len(rows)} 行 / {n_uniq_src} 个 sample")
A(f"删（非队列）      −{len(dropped)} 行（表2 排除名单 29 株中进过 OMArk 的那 10 株）")
A(f"保留              {len(out_rows)} 行  ← 当前 202 株标准")
A("```")
A("")
A("其中保留的 202 行：**192 株**取 `run=original` 行，**10 株**取 `run=filtered` 行"
  "（OMArk 污染过滤后重跑，与 BUSCO 表换 after 值的是同一批）。")
A("")
A(f"> 表2 排除名单是 29 株，为什么这里只删 10 株？——另外 19 株在更早的阶段就被排除了，"
  f"压根没跑 OMArk。对比：BUSCO 213 表删 11 株，比本表多一株 `Rhizopus sp. PT2906`"
  f"（BUSCO D=91.6%，建树前排除，也没跑 OMArk）。")
A("")
A("## 1. 删除的 10 行")
A("")
A("| 原表 sample | 表1编号 | 采纳名 | 问题类别 | 排除原因 | 保留代表 |")
A("|---|---|---|---|---|---|")
for r, e in sorted(dropped, key=lambda x: x[0]['sample']):
    note = ''
    if collapse(r['sample']) == collapse(e['accepted']) and e['orig'] != e['accepted']:
        note = '（原表用采纳名出现）'
    A(f"| `{r['sample']}`{note} | {e['t1_no']} | {e['accepted']} | {e['category']} | "
      f"{str(e['reason']).replace('|', '/')} | {e['keep'] or '—'} |")
A("")
A("> **`Linnemannia_camargensis` 的解释**：它既是现队列 No69 的**采纳名**，又是表2 #209 "
  "（JGI 镜像装配，MIRROR_COPY）的原名。原表里 No69 已由 `Haplosporangium_gracile`（原名）"
  "占位，所以这一条只能判为表2 已排除的镜像记录 —— 看起来像「同一株出现两次」，其实是镜像装配。")
A("")
A("## 2. 重新命名的（第 2 列 现用名）")
A("")
A("第 1 列 `sample` **原样照搬原表**（下划线写法，便于回溯到 OMArk 的 `sum_file` / `run`）；"
  "第 2 列 `现用名` 才是取名结果。")
A("")
A("### 2.1 改名生效：9 株")
A("")
A("| No | 原表 sample | 表3 原名 | 现用名 | 表里看得出变化？ |")
A("|---|---|---|---|---|")
for b in backbone:
    if collapse(b['orig']) != collapse(b['accepted']) and b['no'] not in REVERT_TO_ORIGINAL:
        srcname = claimed[b['no']]['sample']
        curdisp = out_by_no[b['no']]['现用名']
        visible = '是' if srcname != curdisp else '否（原表本就用新名）'
        A(f"| {b['no']} | `{srcname}` | {b['orig']} | `{curdisp}` | {visible} |")
A("")
A("> 其中 No186/187/188/189/192 这 5 株，**原表里用的本来就是采纳名**"
  "（`Actinomortierella_ambigua` 等），所以 `现用名` 与 `sample` 字面相同 —— 不是漏改，"
  "而是原表已经改过了（相对表3 原名仍是改名）。真正在表里字面发生变化的是 "
  "**No12、No57、No58、No69 这 4 株**。")
A("")
A("### 2.2 保留原名：6 株（按你的要求）")
A("")
A("采纳名与表内另一株重名，若采用会把两株并成一株 → 保留原名，因此两列字面相同：")
A("")
for b in backbone:
    if b['no'] in REVERT_TO_ORIGINAL:
        twin = [x for x in backbone if x['accepted'] == b['accepted'] and x['no'] != b['no']]
        who = '、'.join(f"No{t['no']} {t['orig']}" for t in twin)
        A(f"- No{b['no']} **`{b['orig'].replace(' ', '_')}`**（采纳名 {b['accepted']} 与 {who} 重名）")
A("")
A("### 2.3 其余 187 株")
A("")
A("原名与采纳名本就相同，两列一致。**202 个 sample 与 202 个现用名都互不重复。**")
A("")
A("### 2.4 一个写法提醒")
A("")
A("原表的 `sample` 有 6 处与表3 原名的标点写法不同，**只是分隔符差异、不是改名**，"
  "例如 `Absidia_sp.YW-2025a` ↔ 表3 `Absidia sp. YW-2025a`、`Mortierella_sp._KOD1030` ↔ "
  "`Mortierella sp. KOD_1030`。本表 `现用名` 一律沿用原表写法，以免出现「假改名」。"
  "（注：`Mortierella_sp._KOD1030` 与表3 的去分隔符形式一致，其余同理。）")
A("")
A("## 3. 列结构")
A("")
A(f"- 第 1 列 `sample`（原名）、第 2 列 `现用名`（新增）在表头，其余 {len(ORIG_COLS)-1} 列"
  f"（`run, is_final, lineage, conserved_HOGs, proteins, complete_n/pct, S/D/M/A/I/C/U …`）"
  f"原样照搬，未做任何换算。")
A(f"- 共 {len(OUT_COLS)} 列 × {len(out_rows)} 行。")
A("- xlsx 里 `run=filtered` 的 10 行标黄。")
A("")
A("## 4. 校验")
A("")
A(f"- 行数 202；`sample` 唯一 202；`现用名` 唯一 202 ✓")
A(f"- 除新增列外，50 列数值逐格等于原表对应行 ✓（见下方核对）")
A("- 与 BUSCO 202 表的 `sample` 列可以一一对齐（同一套原名）。")
open(REPORT, 'w', encoding='utf-8').write("\n".join(L) + "\n")
print(f"[写出] {os.path.basename(REPORT)}")

# ================================================================== 7. 逐格校验
bad = []
for b in backbone:
    src = claimed[b['no']]
    o = out_by_no[b['no']]
    for c in ORIG_COLS:
        if o[c] != src[c]:
            bad.append((b['no'], c, src[c], o[c]))
assert not bad, bad[:10]
print(f"[校验] 除新增「现用名」列外，其余 {len(ORIG_COLS)-1} 列 × 202 行逐格等于原表 ✓")
print(f"[校验] 删除 {len(dropped)} 行，全部命中表2 排除名单_29 ✓")
print("\n===== 删除清单 =====")
for r, e in sorted(dropped, key=lambda x: x[0]['sample']):
    print(f"   - {r['sample']:<34} [{e['category']}] {str(e['reason'])[:52]}")
print(f"\n===== 改名生效 {n_diff} 株 =====")
for b in backbone:
    if cur_of[b['no']] != b['orig']:
        print(f"   - No{b['no']:>3}  {b['orig']}  ->  {cur_of[b['no']]}")
print("\n===== 保留原名 6 株 =====")
for b in backbone:
    if b['no'] in REVERT_TO_ORIGINAL:
        print(f"   - No{b['no']:>3}  {b['orig']}（采纳名 {b['accepted']} 重名）")
print("\n完成。")
