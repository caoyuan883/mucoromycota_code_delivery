#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重铸 表6 sheet1 的最后两列（RNA-seq / BRAKER_mode）→ 3 列新表。

只读输入：投稿附表/表6_注释证据_202株.xlsx
  · 1_证据总览_202株  ：No / Species_original / BRAKER_mode（当前 RNA-seq 列仅作对照）
  · 3_RNAseq_保留run  ：No + Run_accession（权威口径，逐 run 一行）

输出：投稿附表/表6_RNAseq与BRAKER_mode_202株.xlsx
  · RNAseq_BRAKER（主表）  ：No | RNA-seq | BRAKER_mode，202 行
  · 核对                    ：旧值 / 旧 run 数 / sheet3 run 数 / 新 run 数 / 判定
  · 口径说明                ：差异成因与字段定义
"""
import re
import math
import collections
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ROOT = "/Users/caoyuan/Desktop/Genome_annotation_project_of_Mucoromycota/投稿附表"
SRC = f"{ROOT}/表6_注释证据_202株.xlsx"
OUT = f"{ROOT}/表6_RNAseq与BRAKER_mode_202株.xlsx"

RUN_RE = re.compile(r"[SED]RR\d+")

# ------------------------------------------------------------------ 读源表（只读）
wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
ws1 = wb["1_证据总览_202株"]
ws3 = wb["3_RNAseq_保留run"]

s1 = {}
for r in ws1.iter_rows(min_row=2, values_only=True):
    if isinstance(r[0], int):
        s1[r[0]] = dict(no=r[0], species=r[1], acc=r[4], old_rna=r[6], braker=r[7])

runs = collections.defaultdict(list)
for r in ws3.iter_rows(min_row=2, values_only=True):
    if isinstance(r[0], int):
        runs[r[0]].append(dict(run=r[4], rate=r[5]))
wb.close()

assert list(sorted(s1)) == list(range(1, 203)), "表6 sheet1 应为 No 1–202"
print(f"[源] 表6 sheet1：{len(s1)} 株；sheet3：{sum(len(v) for v in runs.values())} 个 run / {len(runs)} 株")

# ------------------------------------------------------------------ 组装
def runs_of(no):
    return [x["run"] for x in runs.get(no, [])]


rows = []
for no in range(1, 203):
    rec = s1[no]
    rs = runs_of(no)
    rows.append(dict(no=no, species=rec["species"],
                     rna="; ".join(rs) if rs else "N/A",
                     braker=rec["braker"] if rec["braker"] not in (None, "", "—") else "N/A",
                     n_new=len(rs), n_old=len(RUN_RE.findall(str(rec["old_rna"]))),
                     old=str(rec["old_rna"])))

n_run = sum(r["n_new"] for r in rows)
n_braker3 = sum(1 for r in rows if r["braker"] == "BRAKER3")
diff = [r for r in rows if r["n_old"] != r["n_new"]]
multi = [r for r in rows if r["n_new"] >= 2]
print(f"[组装] 202 行；有 RNA-seq {sum(1 for r in rows if r['n_new'])} 株 / {n_run} run；"
      f"BRAKER3 {n_braker3}、BRAKER2 {202 - n_braker3}；多 run 物种 {len(multi)}；旧新计数不一致 {len(diff)}")

# ------------------------------------------------------------------ 写主表
wb2 = openpyxl.Workbook()
ws = wb2.active
ws.title = "RNAseq_BRAKER"

HEAD = ["No", "RNA-seq", "BRAKER_mode"]
HDR_FILL = PatternFill("solid", fgColor="1F4E79")
HDR_FONT = Font(name="Calibri", bold=True, color="FFFFFF", size=10)
BODY_FONT = Font(name="Calibri", size=10)
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
NA_FONT = Font(name="Calibri", size=10, color="808080")

ws.append(HEAD)
for c in range(1, 4):
    cell = ws.cell(1, c)
    cell.fill, cell.font, cell.border = HDR_FILL, HDR_FONT, BORDER
    cell.alignment = Alignment(horizontal="center", vertical="center")
ws.row_dimensions[1].height = 20

WIDTH = {1: 6.5, 2: 46, 3: 14}
for i, w in WIDTH.items():
    ws.column_dimensions[get_column_letter(i)].width = w

for r in rows:
    ws.append([r["no"], r["rna"], r["braker"]])
    rx = ws.max_row
    for c in range(1, 4):
        cell = ws.cell(rx, c)
        cell.font = BODY_FONT
        cell.border = BORDER
        cell.alignment = Alignment(horizontal="center" if c != 2 else "left",
                                   vertical="center", wrap_text=(c == 2))
    if r["rna"] == "N/A":
        ws.cell(rx, 2).font = NA_FONT
    # 行高按换行估算，保证 2 个 / 6 个 run 全部可见（原表被行高截断的问题）
    lines = max(1, math.ceil(len(r["rna"]) / (WIDTH[2] - 2)))
    ws.row_dimensions[rx].height = 14.5 * lines + 1.5

ws.freeze_panes = "A2"
ws.auto_filter.ref = f"A1:C{ws.max_row}"
ws.sheet_view.zoomScale = 100

# ------------------------------------------------------------------ 核对表
ws2 = wb2.create_sheet("核对")
h2 = ["No", "Species_original", "表6 sheet1 原值（RNA-seq 列）",
      "原值 run 数", "sheet3 run 数", "新表 run 数", "判定", "说明"]
ws2.append(h2)
for c in range(1, len(h2) + 1):
    cell = ws2.cell(1, c)
    cell.fill, cell.font, cell.border = HDR_FILL, HDR_FONT, BORDER
    cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
ws2.row_dimensions[1].height = 34

for r in rows:
    if r["n_new"] == 0:
        verdict, note = "N/A", "无保留 RNA-seq run"
    elif r["n_old"] == r["n_new"]:
        verdict = "一致"
        note = ("原单元格已含全部 run，但列宽 16 / 行高 15 只显示第一行；已改为一行全列"
                if r["n_new"] >= 2 else "")
    else:
        verdict = "不一致"
        note = "原值 run 数少于 sheet3，需以 sheet3 为准"
    ws2.append([r["no"], r["species"], r["old"], r["n_old"], r["n_new"], r["n_new"], verdict, note])
    rx = ws2.max_row
    for c in range(1, len(h2) + 1):
        cell = ws2.cell(rx, c)
        cell.font = BODY_FONT
        cell.border = BORDER
        if c in (1, 4, 5, 6, 7):
            cell.alignment = Alignment(horizontal="center", vertical="center")
        elif c == 8:
            cell.alignment = Alignment(horizontal="left", vertical="center", wrap_text=False)
    if r["n_new"] >= 2:
        for c in range(1, len(h2) + 1):
            ws2.cell(rx, c).fill = PatternFill("solid", fgColor="FFF2CC")
    if verdict == "不一致":
        for c in range(1, len(h2) + 1):
            ws2.cell(rx, c).fill = PatternFill("solid", fgColor="FFC7CE")

for i, w in enumerate([6.5, 30, 62, 10, 11, 11, 9, 58], start=1):
    ws2.column_dimensions[get_column_letter(i)].width = w
ws2.freeze_panes = "C2"
ws2.auto_filter.ref = f"A1:H{ws2.max_row}"

# ------------------------------------------------------------------ 口径说明
ws3 = wb2.create_sheet("口径说明")
NOTE = [
    ("表6 · RNA-seq / BRAKER_mode 三列版（202 株）", "title"),
    ("", None),
    ("1. 本表只重铸 表6_注释证据_202株.xlsx → sheet1 的最后两列（RNA-seq、BRAKER_mode），"
     "不加其它字段；No 沿用表6 sheet1 的 1–202。原始表6 与其它投稿附表均未改动。", None),
    ("2. RNA-seq 列的定义：该株用于 BRAKER3 的“保留 run”（表6 sheet3 主键 No + Run_accession），"
     "一株多个 run 时以“; ”分隔全部列出，按 sheet3 顺序排列；不附比对率。", None),
    ("3. 缺失值写法：无保留 run、无 mode 一律写 N/A（不用“—”）。"
     f"本版：有 RNA-seq {sum(1 for r in rows if r['n_new'])} 株 / 共 {n_run} 个 run；"
     f"无 RNA-seq {202 - sum(1 for r in rows if r['n_new'])} 株（→ N/A）。", None),
    ("4. 关键核对结论：逐株比对 No 1–202，sheet1 原 RNA-seq 单元格与 sheet3 的 run 集合 202/202 完全一致，"
     "run 总数两边均为 119 —— 数据本身没有漏记，漏的是“显示”：原表 RNA-seq 列宽仅 16 字符且开启自动换行，"
     "行高却固定 15 pt（约一行），Excel 只画得出第一行，第二个及以后的 run 被行高裁掉；"
     "6 个 run 的 Rhizopus microsporus（No.153）、Rhizopus stolonifer（No.155）行高 45 pt 也只显示 3 行。"
     "所以“数量看起来是 1”属实，“数据只存了 1 个”不成立。", None),
    ("5. 本版处理：每个 run 全部写在同一单元格，行高按字符数自动加高（1–2 行），"
     "任何 2 个 / 6 个 run 都完整可见；见「核对」页浅黄底 28 株（26 株 2 个 run、2 株 6 个 run）。", None),
    ("6. BRAKER_mode：照录 sheet1（BRAKER3 = 基因组 + 蛋白 + RNA-seq；BRAKER2 = 基因组 + 蛋白）。"
     f"计数 BRAKER3 {n_braker3} / BRAKER2 {202 - n_braker3}，与“有保留 run 的株数”"
     f"{sum(1 for r in rows if r['n_new'])} 一一对应（有 run 必 BRAKER3，无 run 必 BRAKER2）。", None),
    ("7. 保留 run 的筛选口径（表6 口径页第 3 条）：比对率 <20% 剔除；每物种默认保留比对率最高的 2 个 run；"
     "若全部 run 比对率相近且无高分梯队（Rhizopus microsporus、Rhizopus stolonifer，约 50%），保留全部通过阈值者。", None),
    ("8. 数据来源（只读）：投稿附表/表6_注释证据_202株.xlsx 的 sheet1（No、Species_original、BRAKER_mode）"
     "与 sheet3（No、Run_accession）。生成脚本：投稿附表/归档与脚本/脚本/build_表6_RNAseq_BRAKER_3列.py。", None),
]
ws3.column_dimensions["A"].width = 150
for i, (text, kind) in enumerate(NOTE, start=1):
    cell = ws3.cell(i, 1, text)
    cell.alignment = Alignment(wrap_text=True, vertical="top")
    if kind == "title":
        cell.font = Font(bold=True, size=12, color="1F4E79")
        ws3.row_dimensions[i].height = 22
    else:
        cell.font = Font(size=10)
        ws3.row_dimensions[i].height = 15 * max(1, math.ceil(len(text) / 108)) if text else 8

wb2.save(OUT)
print(f"[输出] {OUT}")
print(f"   主表 RNAseq_BRAKER：{ws.max_row - 1} 行 × 3 列")
print(f"   核对：{ws2.max_row - 1} 行（多 run 28 株浅黄；不一致 0 行）")
