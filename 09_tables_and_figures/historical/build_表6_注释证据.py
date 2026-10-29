#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""表6 注释证据表 v3（命名对齐表3/表5 + Sheet3 主键）装配与输出（只读源文件，输出新文件）
v3 修改：Sheet1 增 Species_original 主名列（原名），原 Species 列改名 Species_current_IF；
        Sheet2 表A 同步两列命名；Sheet3 前置 No + Assembly_accession（主键 No + Run_accession），
        移除非队列 Amylomyces rouxii / Podila humilis 两行。
"""
import openpyxl, csv, re
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from collections import defaultdict

ROOT = "/Users/caoyuan/Desktop/Genome_annotation_project_of_Mucoromycota"
OUT = f"{ROOT}/投稿附表/表6_注释证据_202株.xlsx"

def collapse(name):
    return re.sub(r'[\s_.\-]+', '', str(name)).lower()

def italic_rule(name):
    s = str(name)
    if 'sp.' in s or ' ' not in s:
        return False
    parts = s.split(' ')
    return len(parts) == 2 and all('.' not in p for p in parts)

# ================================================================ 表3 骨干
wb = openpyxl.load_workbook(f"{ROOT}/投稿附表/表3_研究总表.xlsx", read_only=True)
ws = wb["研究菌株_202"]
backbone = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == '编号':
        continue
    backbone.append(dict(no=r[0], orig=r[1], strain=r[2], accepted=r[3],
                         acc=r[4], db=r[5], biosample=r[6]))
wb.close()

# ================================================================ 蛋白证据评估（219 样本）
wb = openpyxl.load_workbook(f"{ROOT}/6_Protein_evidence/Selection_of_protein_evidence.xlsx", read_only=True)
ws = wb["Sheet1"]
selection = []
for r in ws.iter_rows(values_only=True):
    if r[1] is None or r[1] == 'sample':
        continue
    selection.append(dict(sample=r[1], subphylum=r[2], cls=r[3], order=r[4],
                          family=r[5], genus=r[6], species=r[7],
                          mucoro=r[8], fungi=r[9], lib=r[10], idx=r[11]))
wb.close()

# ================================================================ OrthoDB157（133 属种，出现次数合计 157）
wb = openpyxl.load_workbook(f"{ROOT}/6_Protein_evidence/14_tables_OrthoDBv12.2_Mucoromycota157_ProteinEvidence_Taxonomy.xlsx", read_only=True)
ws = wb["Sheet1"]
ortho157 = []
for r in ws.iter_rows(values_only=True):
    if r[0] is None or r[0] == 'Index':
        continue
    ortho157.append(dict(sample=r[1], subphylum=r[2], cls=r[3], order=r[4],
                         family=r[5], genus=r[6], species=r[7], occ=r[8]))
wb.close()

# ================================================================ 模式映射（211）
braker3, braker2 = {}, set()
for line in open(f"{ROOT}/7_RNA-seq_data_preparation/species_mode_mapping.txt"):
    line = line.strip()
    if not line:
        continue
    m = re.match(r'BRAKER3: (\S+) -> (.*)$', line)
    if m:
        sp = m.group(1)
        runs = sorted(set(re.search(r'__((?:SRR|ERR|DRR)\w+)\.sorted\.bam$', b).group(1)
                          for b in m.group(2).split()))
        braker3[sp] = runs
    else:
        m2 = re.match(r'BRAKER2: (\S+) \(无RNA-seq\)', line)
        if m2:
            braker2.add(m2.group(1))
mapping_tags = set(braker3) | braker2

# ================================================================ QC 21 修订表
wb = openpyxl.load_workbook(f"{ROOT}/7_RNA-seq_data_preparation/21_tables_revised_rnaseq_qc_summary.xlsx", read_only=True)
ws = wb["工作表 1 - rnaseq_qc_summary"]
qc_main = {}
for r in ws.iter_rows(values_only=True):
    if r[0] and r[0] != 'Sample':
        qc_main[r[0]] = r[1:9]
ws = wb["Sheet1"]
qc_excl = {}
for r in ws.iter_rows(values_only=True):
    if r[0] and r[0] != 'Sample' and not str(r[0]).startswith('Note'):
        qc_excl[r[0]] = (r[1], r[9])
wb.close()

rename_map = {
    'Mortierella_capitata': 'Actinomortierella_capitata',
    'Mortierella_wolfii': 'Actinomortierella_wolfii',
    'Mortierella_ambigua': 'Actinomortierella_ambigua',
    'Mortierella_elongata': 'Linnemannia_elongata',
    'Mortierella_gamsii': 'Linnemannia_gamsii',
    'Mortierella_humilis': 'Podila_humilis',
    'Mortierella_minutissima': 'Podila_minutissima',
    'Mortierella_multidivaricata': 'Gamsiella_multidivaricata',
    'Mortierella_verticillata': 'Podila_verticillata',
}
def canon(sp):
    return rename_map.get(sp, sp)

run_rows = {}
for key, vals in qc_main.items():
    sp0, run = key.split('__', 1)
    sp = canon(sp0)
    run_rows[(sp, run)] = dict(species=sp, run=run,
        rate=vals[0], total=vals[1], mapped=vals[2], paired=vals[3],
        duprate=vals[4], afterfilter=vals[5], insert=vals[6], bamsize=vals[7])
for key, (rate, cat) in qc_excl.items():
    sp0, run = key.split('__', 1)
    sp = canon(sp0)
    k = (sp, run)
    if k not in run_rows:
        run_rows[k] = dict(species=sp, run=run, rate=rate,
            total=None, mapped=None, paired=None, duprate=None,
            afterfilter=None, insert=None, bamsize=None)

# ================================================================ QC 账本
qc_pipeline = {}
for r in csv.DictReader(open(f"{ROOT}/证据链归档/09_蛋白组BUSCO评估/结果_主表/final_protein_qc_master_table.csv")):
    qc_pipeline[r['sample']] = r['annotation_pipeline']

# ================================================================ 别名与匹配
alias_map = {
    'Absidia_sp._YW-2025a': 'Absidia_sp.YW-2025a',
    'Absidia_sp._YW-2025b': 'Absidia_sp.YW-2025b',
    'Entrophospora_sp._SA101': 'Entrophospora_sp.SA101',
    'Mortierella_ambigua': 'Actinomortierella_ambigua',
    'Mortierella_capitata': 'Actinomortierella_capitata',
    'Mortierella_wolfii': 'Actinomortierella_wolfii',
    'Mortierella_minutissima': 'Podila_minutissima',
    'Mortierella_multidivaricata': 'Gamsiella_multidivaricata',
    'Mortierella_sp._BCC40632': 'Mortierella_sp.BCC40632',
    'Mortierella_sp._KOD_1030': 'Mortierella_sp._KOD1030',
    'Rhizophagus_sp._LPA8_1_AMF': 'Rhizophagus_sp.LPA8_1_AMF',
    'Rhizopus_sp._FBL578': 'Rhizopus_sp.FBL578',
}
def tag_of(name):
    return str(name).replace(' ', '_')

def find_mapping_tag(orig, accepted):
    cands = [tag_of(orig), tag_of(accepted),
             alias_map.get(tag_of(orig)), alias_map.get(tag_of(accepted))]
    for c in cands:
        if c and c in mapping_tags:
            return c
    ccol = collapse(orig)
    for t in mapping_tags:
        if collapse(t) == ccol:
            return t
    return None

sel_collapse = {collapse(s['sample']): s for s in selection}
ortho_collapse = {collapse(s['sample']): s for s in ortho157}
qc_collapse = {collapse(k): v for k, v in qc_pipeline.items()}

genus_tax = {
    'Mortierella': ('Mortierellomycotina', 'Mortierellomycetes'),
    'Cryptocalciella': ('Mortierellomycotina', 'Mortierellomycetes'),
    'Mucor': ('Mucoromycotina', 'Mucoromycetes'),
    'Mycotypha': ('Mucoromycotina', 'Mucoromycetes'),
    'Syzygites': ('Mucoromycotina', 'Mucoromycetes'),
    'Umbelopsis': ('Mucoromycotina', 'Umbelopsidomycetes'),
    'Rhizopus': ('Mucoromycotina', 'Mucoromycetes'),
}

master = []
for b in backbone:
    mtag = find_mapping_tag(b['orig'], b['accepted'])
    cc = collapse(b['orig'])
    sel = sel_collapse.get(cc) or sel_collapse.get(collapse(b['accepted'])) or \
          (sel_collapse.get(collapse(mtag)) if mtag else None)
    ortho = ortho_collapse.get(cc) or (ortho_collapse.get(collapse(mtag)) if mtag else None)
    qm = qc_collapse.get(cc) or qc_collapse.get(collapse(b['accepted'])) or \
         (qc_collapse.get(collapse(mtag)) if mtag else None)
    row = dict(backbone=b, mtag=mtag, sel=sel, ortho=ortho, qc_mode=qm)
    row['rna_seq'] = 'No' if (mtag is None or mtag in braker2) else 'Yes'
    if mtag and mtag in braker3:
        row['runs'] = sorted(braker3[mtag],
            key=lambda r: -(run_rows.get((mtag, r), {}).get('rate') or -1))
    else:
        row['runs'] = []
    row['cand'] = [(k[1], run_rows[k]) for k in run_rows if k[0] == mtag] if mtag else []
    master.append(row)

for r in master:
    if r['sel'] and r['sel']['subphylum']:
        continue
    genus = str(r['backbone']['orig']).split(' ')[0]
    r['tax_filled'] = genus_tax.get(genus)

# 采纳名相同的组（用于 Note 标注）
acc_groups = defaultdict(list)
for r in master:
    acc_groups[r['backbone']['accepted']].append(r['backbone']['no'])

# ================================================================ 统计
n_yes = sum(1 for r in master if r['rna_seq'] == 'Yes')
n_no = sum(1 for r in master if r['rna_seq'] == 'No')
n_covered = sum(1 for r in master if r['ortho'])
n_uncovered = len(master) - n_covered
n_retained = sum(len(r['runs']) for r in master)
print(f"队列: 202 | BRAKER3={n_yes} BRAKER2={n_no} | 保留run={n_retained} | "
      f"OrthoDB157覆盖={n_covered} 未覆盖={n_uncovered}")
print("同名采纳名组:", {k: v for k, v in acc_groups.items() if len(v) > 1})
print("未映射:", [(r['backbone']['no'], r['backbone']['orig'], r['qc_mode']) for r in master if r['mtag'] is None])

# ================================================================ 写 xlsx
wb = openpyxl.Workbook()
hdr_font = Font(bold=True, color="FFFFFF", size=10)
hdr_fill = PatternFill("solid", fgColor="1F4E79")
hdr_align = Alignment(horizontal="center", vertical="center", wrap_text=True)
thin = Side(style="thin", color="BFBFBF")
border = Border(left=thin, right=thin, top=thin, bottom=thin)
body_font = Font(size=10)
wrap = Alignment(vertical="top", wrap_text=True)
note_font = Font(size=10, color="404040")
bold_note = Font(size=10, bold=True)

def style_sheet(ws, ncols, widths, freeze="A2"):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=1, column=c)
        cell.font = hdr_font
        cell.fill = hdr_fill
        cell.alignment = hdr_align
        cell.border = border
        ws.column_dimensions[get_column_letter(c)].width = widths[c - 1]
    ws.freeze_panes = freeze
    ws.auto_filter.ref = f"A1:{get_column_letter(ncols)}{ws.max_row}"
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, max_col=ncols):
        for cell in row:
            if not (cell.font and cell.font.i) and not cell.font.bold:
                cell.font = body_font
            cell.alignment = wrap
            cell.border = border

def set_species_cell(cell, name):
    cell.value = name
    if italic_rule(name):
        cell.font = Font(size=10, italic=True)

# ================= Sheet 1 证据总览（命名对齐表3/表5）
ws1 = wb.active
ws1.title = "1_证据总览_202株"
cols1 = ["No", "Species_original", "Species_current_IF", "Assembly_accession",
         "Protein", "RNA-seq", "BRAKER_mode"]
ws1.append(cols1)
for r in master:
    b, mtag = r['backbone'], r['mtag']
    if r['runs']:
        rna = "; ".join("%s (%.2f%%)" % (a, float(run_rows[(mtag, a)]['rate']))
                        for a in r['runs'])
    else:
        rna = "—"
    mode = "BRAKER3" if r['rna_seq'] == 'Yes' else "BRAKER2"
    ws1.append([b['no'], b['orig'], b['accepted'], b['acc'], "+", rna, mode])
    set_species_cell(ws1.cell(row=ws1.max_row, column=2), b['orig'])
    set_species_cell(ws1.cell(row=ws1.max_row, column=3), b['accepted'])
style_sheet(ws1, len(cols1), [5, 24, 24, 18, 8, 52, 11])

# ---- Sheet1 底部说明
first_note = ws1.max_row + 2
notes1 = [
    ("表6 · 注释证据表（说明）", True),
    ("", False),
    ("• 物种命名：Species_original 为原名（主名），Species_current_IF 为采纳名——与表3、表5 及物种树策略一致"
     "（统一用原名 + 附对照表）。采纳名相同的 5 组 11 行（Cunninghamella elegans、Linnemannia hyalina、"
     "Mucor circinelloides、Rhizopus arrhizus、Rhizopus microsporus）按 No / Assembly_accession 区分。", False),
    ("• 蛋白证据（“+”）：202 株全部使用 OrthoDB v12.2 的 Fungi 分区蛋白集（2,692 个真菌蛋白质文件；"
     "BRAKER2/3 均经 --prot_seq 使用清理后的 Fungi.clean.fa）作为统一的蛋白质证据。", False),
    ("• RNA-seq：括号内为 HISAT2 overall alignment rate（%）。每物种默认保留比对率最高的 2 个 run；"
     "当该物种全部候选 run 比对率相近、无明显高分梯队时（Rhizopus microsporus 6 个 run 43.93–62.42%、"
     "Rhizopus stolonifer 6 个 run 47.87–66.33%），保留全部通过阈值的 run；比对率 <20% 的 run 一律剔除"
     "（Geosiphon pyriformis、Ambispora leptoticha、Paraglomus brasilianum 的候选 run 因此全部被剔除，无 RNA-seq 可用，记为 “—”）。", False),
    ("• BRAKER 模式：BRAKER3 = 基因组 + 蛋白证据 + RNA-seq；BRAKER2 = 基因组 + 蛋白证据（仅蛋白）。", False),
    ("• 合计 202 株：BRAKER3（蛋白 + RNA-seq）83 株；BRAKER2（仅蛋白）119 株。队列内保留 RNA-seq run 共 119 个"
     "（见 Sheet3，主键为 No + Run_accession）。", False),
    ("• 本表结构参考 Nenasheva et al. 2025, Scientific Data 12:985 的 Table 1（每物种的基因组与转录组数据来源）。", False),
]
for text, bold in notes1:
    ws1.cell(row=first_note, column=1, value=text)
    c = ws1.cell(row=first_note, column=1)
    c.font = bold_note if bold else note_font
    c.alignment = Alignment(wrap_text=True, vertical="top")
    ws1.merge_cells(start_row=first_note, start_column=1, end_row=first_note, end_column=7)
    first_note += 1

# ================= Sheet 2 蛋白证据
ws2 = wb.create_sheet("2_蛋白证据_OrthoDB")
for colw, wdt in zip("ABCDEFGHIJ", [8, 26, 24, 18, 16, 16, 30, 20, 12, 10]):
    ws2.column_dimensions[colw].width = wdt
r0 = 1
secA = [
    ("蛋白证据说明（研究记录：6_Protein_evidence/（OrhoDB V12.2 公共蛋白质文件研究）*.md）", True),
    ("• OrthoDB v12.2 层级蛋白库规模：Fungi 2,692 个蛋白质文件；Mucoromycota 157 个蛋白质文件（对应 133 个不同属种）；Mucorales 57 个蛋白质文件。", False),
    ("• 157 文件的亚门分布：Glomeromycotina 26 / Mortierellomycotina 48 / Mucoromycotina 59（合计 133 属种，同属种出现次数合计 157；"
     "重复较多者如 Rhizophagus irregularis ×7、Mortierella alpina ×5、Rhizopus arrhizus ×3）。", False),
    ("• 覆盖检查：202 队列中仅 %d 株在 Mucoromycota 专属库中有对应属种（%d 株未覆盖，见下表 F 列，可筛选 “No”），" % (n_covered, n_uncovered) +
     "因此最终选用覆盖更广的 Fungi 级作为统一蛋白证据。", False),
    ("• 论文表述要点：为确保跨基因组统一注释时各蛋白家族尽可能得到代表，采用 OrthoDB v12.2 的 Fungi 分区蛋白集（由 BRAKER 维护团队提供、可直接用于 BRAKER2/3）；"
     "虽然 Mucoromycota 专属蛋白集在系统发育上更近，但考虑到统一流程批量注释，优先选择覆盖更广的 Fungi 数据集，以降低近缘参考不足导致的 hints 缺失。", False),
    ("", False),
    ("A. 202 队列覆盖检查（OrthoDB v12.2 Mucoromycota 157 库）", True),
]
for text, bold in secA:
    ws2.cell(row=r0, column=1, value=text)
    c = ws2.cell(row=r0, column=1)
    c.font = bold_note if bold else note_font
    c.alignment = Alignment(wrap_text=True, vertical="top")
    ws2.merge_cells(start_row=r0, start_column=1, end_row=r0, end_column=7)
    r0 += 1

hdrB = ["No", "Species_original", "Species_current_IF", "Subphylum", "Class",
        "OrthoDB_Mucoromycota157", "Note"]
hdr_row = r0
for j, h in enumerate(hdrB, 1):
    cell = ws2.cell(row=r0, column=j, value=h)
    cell.font = hdr_font
    cell.fill = hdr_fill
    cell.alignment = hdr_align
    cell.border = border
r0 += 1
tbl_start = r0
for r in master:
    b, ortho = r['backbone'], r['ortho']
    sub = (r['sel'] or {}).get('subphylum') or (r.get('tax_filled') or ('—',))[0]
    cls = (r['sel'] or {}).get('cls') or (r.get('tax_filled') or ('—', '—'))[1] or "—"
    ortho_s = ("Yes (×%s)" % ortho['occ']) if ortho else "No"
    note = ""
    siblings = [n for n in acc_groups.get(b['accepted'], []) if n != b['no']]
    if siblings:
        note = "采纳名与 No.%s 相同（按 Accession 区分）" % "/".join(str(n) for n in sorted(siblings))
    ws2.cell(row=r0, column=1, value=b['no'])
    ws2.cell(row=r0, column=2, value=b['orig'])
    ws2.cell(row=r0, column=3, value=b['accepted'])
    ws2.cell(row=r0, column=4, value=sub)
    ws2.cell(row=r0, column=5, value=cls)
    ws2.cell(row=r0, column=6, value=ortho_s)
    ws2.cell(row=r0, column=7, value=note if note else "—")
    set_species_cell(ws2.cell(row=r0, column=2), b['orig'])
    set_species_cell(ws2.cell(row=r0, column=3), b['accepted'])
    r0 += 1
tbl_end = r0 - 1
for row in ws2.iter_rows(min_row=tbl_start, max_row=tbl_end, max_col=7):
    for cell in row:
        cell.font = body_font
        cell.alignment = wrap
        cell.border = border
ws2.auto_filter.ref = f"A{hdr_row}:G{tbl_end}"

r0 += 1
secC = "B. OrthoDB v12.2 Mucoromycota 157 库 · 133 属种明细（同属种出现次数合计 157；In_cohort 为该属种是否命中 202 队列）"
ws2.cell(row=r0, column=1, value=secC)
c = ws2.cell(row=r0, column=1)
c.font = bold_note
c.alignment = Alignment(wrap_text=True, vertical="top")
ws2.merge_cells(start_row=r0, start_column=1, end_row=r0, end_column=10)
r0 += 1

hdrC = ["Index", "Sample", "Subphylum", "Class", "Order", "Family", "Genus", "Species",
        "Number_of_occurrences", "In_cohort_202"]
hdr_rowC = r0
for j, h in enumerate(hdrC, 1):
    cell = ws2.cell(row=r0, column=j, value=h)
    cell.font = hdr_font
    cell.fill = hdr_fill
    cell.alignment = hdr_align
    cell.border = border
r0 += 1
no_by_collapse = {}
for r in master:
    no_by_collapse.setdefault(collapse(r['backbone']['orig']), r['backbone']['no'])
    no_by_collapse.setdefault(collapse(r['backbone']['accepted']), r['backbone']['no'])
tblC_start = r0
for i, s in enumerate(ortho157, 1):
    incoh = no_by_collapse.get(collapse(s['sample']))
    ws2.cell(row=r0, column=1, value=i)
    ws2.cell(row=r0, column=2, value=s['sample'])
    ws2.cell(row=r0, column=3, value=s['subphylum'])
    ws2.cell(row=r0, column=4, value=s['cls'])
    ws2.cell(row=r0, column=5, value=s['order'])
    ws2.cell(row=r0, column=6, value=s['family'])
    ws2.cell(row=r0, column=7, value=s['genus'])
    ws2.cell(row=r0, column=8, value=s['species'])
    ws2.cell(row=r0, column=9, value=s['occ'])
    ws2.cell(row=r0, column=10, value=("Y(No.%s)" % incoh) if incoh else "N")
    set_species_cell(ws2.cell(row=r0, column=8), s['species'])
    r0 += 1
tblC_end = r0 - 1
for row in ws2.iter_rows(min_row=tblC_start, max_row=tblC_end, max_col=10):
    for cell in row:
        cell.font = body_font
        cell.alignment = wrap
        cell.border = border
ws2.auto_filter.ref = f"A{hdr_rowC}:J{tblC_end}"

# ================= Sheet 3 RNA-seq 保留 run（仅队列；主键 No + Run_accession）
ws3 = wb.create_sheet("3_RNAseq_保留run")
cols3 = ["No", "Assembly_accession", "Species_original", "Run_accession",
         "AlignmentRate_pct", "Total_reads", "Mapped_reads", "Properly_paired",
         "Reads_after_filter", "Insert_size_peak", "BAM_size_MB"]
ws3.append(cols3)

def num(v):
    if v is None or v == 'NA' or v == '':
        return 'NA'
    try:
        f = float(v)
        return round(f, 2) if f != int(f) else int(f)
    except (TypeError, ValueError):
        return v

n_run_rows = 0
for r in sorted(master, key=lambda x: x['backbone']['no']):
    b, mtag = r['backbone'], r['mtag']
    if not r['runs']:
        continue
    for run in r['runs']:
        rec = run_rows.get((mtag, run), {})
        ws3.append([b['no'], b['acc'], b['orig'], run,
                    num(rec.get('rate')), num(rec.get('total')),
                    num(rec.get('mapped')), num(rec.get('paired')),
                    num(rec.get('afterfilter')), num(rec.get('insert')),
                    num(rec.get('bamsize'))])
        set_species_cell(ws3.cell(row=ws3.max_row, column=3), b['orig'])
        n_run_rows += 1
style_sheet(ws3, len(cols3), [5, 18, 26, 14, 11, 12, 12, 12, 12, 10, 9])

first_note3 = ws3.max_row + 2
notes3 = [
    ("RNA-seq 保留 run 说明（仅列 202 队列内、用于实验的 run；主键：No + Run_accession）", True),
    ("• 重名物种按 No 区分：如 Cunninghamella elegans 的 2 个 run 属 No.41（原名 Cunninghamella bertholletiae）；"
     "Mucor circinelloides 的 4 个 run 分属 No.114（ERR11901554/55）与 No.127（原名 Mucor velutinosus，SRR28755641/42）；"
     "Rhizopus arrhizus 的 4 个 run 分属 No.150（原名 Rhizopus arrhizus）与 No.152（原名 Rhizopus delemar）。", False),
    ("• 筛选规则（研究记录：7_RNA-seq_data_preparation/（RNA-seq + protein 的 BRAKER3_ETP 路线可行性研究）*.md 最终定稿）："
     "比对率 <20% 的 run 一律剔除；每物种默认保留比对率最高的 2 个 run；"
     "当物种全部候选 run 比对率相近、无明显高分梯队时（Rhizopus microsporus、Rhizopus stolonifer，6 个 run 均约 50%），保留全部通过阈值的 run。", False),
    ("• 队列内合计：119 个 run / 83 株。非队列基因组（Amylomyces rouxii、Podila humilis）的保留 run 不列入本表；"
     "全量含排除的 run 级 QC 见归档版 表6 v1/v2 或 21_tables_revised_rnaseq_qc_summary.xlsx。", False),
]
for text, bold in notes3:
    ws3.cell(row=first_note3, column=1, value=text)
    c = ws3.cell(row=first_note3, column=1)
    c.font = bold_note if bold else note_font
    c.alignment = Alignment(wrap_text=True, vertical="top")
    ws3.merge_cells(start_row=first_note3, start_column=1, end_row=first_note3, end_column=11)
    first_note3 += 1

# ================= Sheet 4 口径与来源
ws4 = wb.create_sheet("4_口径与来源")
ws4.column_dimensions['A'].width = 140
lines = [
    ("表6 · 注释证据表 v3 — 规则口径与数据来源", True),
    ("", False),
    ("1. 骨干与命名：表3 研究菌株_202（编号 1–202）。Species_original（原名）为主名，Species_current_IF 为采纳名，"
     "与表3、表5 及物种树策略一致（统一用原名 + 附对照表）；Sheet3 主键为 No + Run_accession。", False),
    ("", False),
    ("2. 蛋白证据（研究记录：6_Protein_evidence/（OrhoDB V12.2 公共蛋白质文件研究）*.md）：", False),
    ("   a) OrthoDB v12.2 层级库：Fungi 2,692 文件 / Mucoromycota 157 文件（133 属种）/ Mucorales 57 文件。", False),
    ("   b) 157 文件亚门分布 26/48/59（Glo/Mor/Mur）；202 队列中 120 株被覆盖、82 株未覆盖（见 Sheet2 表A）。", False),
    ("   c) 结论：为保证统一流程注释的蛋白家族代表性，最终选用 Fungi 级（2,692 文件）作为全部菌株的蛋白证据；"
     "BRAKER 输入文件为清理后的 Fungi.clean.fa（--prot_seq）。", False),
    ("", False),
    ("3. RNA-seq 证据（研究记录：7_RNA-seq_data_preparation/（RNA-seq + protein 的 BRAKER3_ETP 路线可行性研究）*.md）：", False),
    ("   a) 检索与下载：每物种按 \"[Organism] AND illumina[Platform] AND \"rna seq\"[Strategy]\" 查询 SRA；只保留 PAIRED 文库，"
     "每物种最多取前 6 个 run（仿硅藻论文 \"top six entries\"；原始 1,552 条记录 → sci-grade 211 条/91 物种 → 最终 96 物种/215 run 下载清单）。", False),
    ("   b) 处理：fasterq-dump（SRA Toolkit 3.4.1）→ fastp 0.23.4 → HISAT2 比对 → SAMtools（sort/index/flagstat）。", False),
    ("   c) 筛选（最终定稿规则）：比对率 <20% 剔除；每物种默认保留比对率最高的 2 个 run；"
     "若全部 run 比对率相近且无高分梯队（Rhizopus microsporus、Rhizopus stolonifer，均约 50%），保留全部通过阈值的 run。"
     "（阈值经 40%→30%→20% 讨论后敲定：20% 门槛足以剔除明显异常 run，同时最大限度保住物种的 BRAKER3 资格。）", False),
    ("   d) 保留 run 合并为物种级 merged BAM → 有 BAM 者 BRAKER3（基因组+蛋白+RNA-seq，--fungus --busco_lineage mucoromycota_odb12 --softmasking --gff3），无 BAM 者 BRAKER2（基因组+蛋白）。", False),
    ("", False),
    ("4. 数量口径：", False),
    ("   • mode-mapping 211 个软掩蔽基因组（BRAKER3 85 + BRAKER2 126）；其中 200 个对应 202 队列菌株。", False),
    ("   • 最终 202 队列：BRAKER3 83 株 + BRAKER2 119 株；队列内保留 run 119 个（Sheet3 仅列队列内 run；"
     "非队列基因组 Amylomyces rouxii、Podila humilis 的保留 run 已按要求移除）。", False),
    ("   • 有 SRA 候选但全部 run <20% 被剔除（→BRAKER2）的 3 株：Geosiphon pyriformis（6.65/0.04/0.04%）、"
     "Ambispora leptoticha（19.31%）、Paraglomus brasilianum（11.76%）。", False),
    ("   • No.201 Cryptocalciella humilis（GCA_057980025.1）与 No.202 Rhizopus sp. TYR1（GCA_057662165.1）："
     "mode-mapping（07-14）之后加入，据 09 蛋白 QC 账本均为 BRAKER2。", False),
    ("   • 09 蛋白 QC 账本（213 行）= 127 BRAKER2 + 83 BRAKER3 + 3 UNKNOWN（Actinomortierella ambigua / Benjaminiella poitrasii / "
     "Dichotomocladium robustum）；与 mapping 的差异：Rhizophagus diaphanus、Amylomyces rouxii 多跑 1 次 BRAKER3（不入队列）。", False),
    ("", False),
    ("5. 特殊记录：", False),
    ("   • Oehlia diaphana（No.129）的 SRR916888 原按 Rhizophagus diaphanus 检索，后用于本株（比对率 98.08%）。", False),
    ("   • 分类学改名株的 QC 记录含旧名（Mortierella_*），已合并到新名（Actinomortierella/Gamsiella/Podila/Linnemannia）；"
     "21 修订表排除类别计数：Low alignment Rate 15 / Replicate redundancy 66 / Taxonomic redundancy 12。", False),
    ("   • Endogone sp. FLAS-F59071__SRR10513696 仅见 raw QC tsv（物种不在队列）。", False),
    ("", False),
    ("6. 数据来源（全部只读）：", False),
    ("   • 投稿附表/表3_研究总表.xlsx（研究菌株_202）；表5_软屏蔽QC_202株.xlsx（软掩蔽交付清单）", False),
    ("   • 6_Protein_evidence/：Selection_of_protein_evidence.xlsx、14_tables_OrthoDBv12.2_Mucoromycota157_ProteinEvidence_Taxonomy.xlsx、"
     "mucoromycota_157_species_taxonomy.csv、（OrhoDB V12.2 公共蛋白质文件研究）*.md/*.docx", False),
    ("   • 7_RNA-seq_data_preparation/：21_tables_revised_rnaseq_qc_summary.xlsx、species_mode_mapping.txt（22_tables_species_mode_mapping.xlsx）、"
     "rnaseq_qc_summary.tsv、19_tables_species_level_rnaseq_plan.xlsx、（RNA-seq + protein 的 BRAKER3_ETP 路线可行性研究）*.md/*.docx", False),
    ("   • 证据链归档/09_蛋白组BUSCO评估/结果_主表/final_protein_qc_master_table.csv（QC 账本模式列）", False),
    ("   • braker3_batch_v4/00_config.sh（PROT_DB=Fungi.clean.fa；BRAKER3 参数）", False),
    ("", False),
    ("7. 排版规范（同表3/表5）：属种名斜体；sp. 及其后菌株编号不斜体。", False),
]
for text, bold in lines:
    ws4.append([text])
    c = ws4.cell(row=ws4.max_row, column=1)
    c.font = Font(bold=bold, size=10)
    c.alignment = Alignment(wrap_text=True, vertical="top")

wb.save(OUT)
print("已写出:", OUT)
print("Sheet3 run 行数:", n_run_rows)
