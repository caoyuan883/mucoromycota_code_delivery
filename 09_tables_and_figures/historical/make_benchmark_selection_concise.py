#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成《基准物种遴选 · 简洁版（含 Scientific Data Methods 原文）》。

输入（只读）：
  基准物种20株_名单_v2.csv          —— V2 遴选的 20 株及其设计格
  benchmark_species_evidence.json   —— NCBI 装配与注释元数据（2026-09-14 快照）
输出：
  基准物种遴选_简洁版_含Methods原文.docx

运行：python3 9_Benchmark_reference_annotations/make_benchmark_selection_concise.py
"""
import os
import json
import pandas as pd

OUTDIR = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(OUTDIR, "基准物种20株_名单_v2.csv")
DOCX = os.path.join(OUTDIR, "基准物种遴选_简洁版_含Methods原文.docx")

# 参考注释流程与文献（用于表 1 与英文 Methods；均已核实到一手来源）
REF_PROV = {
    "GCF_026210795.1": ("BRAKER 2.1.5 + funannotate/PASA", "Manley et al. 2023, G3 13:jkad077"),
    "GCF_035048745.1": ("BRAKER1 v1.9 + CodingQuarry + TransDecoder → EVM + funannotate 1.8.13",
                        "Conlan et al. 2024, Microbiol Resour Announc 13:e00062-24"),
    "GCA_015698045.1": ("BRAKER v2.1.2 + AUGUSTUS 3.3.2 + GeneMark-ES 4.32",
                        "Kobayashi et al. 2026, Research Square preprint (under revision)"),
    "GCA_040256725.1": ("JGI Annotation Pipeline", "JGI MycoCosm portal Mucci3"),
    "JGI:Disdec1_v1.0": ("JGI Annotation Pipeline", "JGI MycoCosm portal Disdec1"),
    "GCA_036320915.1": ("JGI Annotation Pipeline", "JGI MycoCosm portal MoeNVP64_1"),
    "GCA_977092135.1": ("EBP-Nor pipeline: miniprot + GALBA + EvidenceModeler",
                        "github.com/ebp-nor/GenomeAnnotation"),
    "GCF_029851405.1": ("BRAKER2, protein-only mode", "Shea et al. 2023, Microbiol Resour Announc 12:e00338-23"),
    "GCA_977110935.1": ("EBP-Nor pipeline: miniprot + GALBA + EvidenceModeler",
                        "github.com/ebp-nor/GenomeAnnotation"),
    "GCA_977110905.1": ("EBP-Nor pipeline: miniprot + GALBA + EvidenceModeler",
                        "github.com/ebp-nor/GenomeAnnotation"),
    "GCA_977110955.1": ("EBP-Nor pipeline: miniprot + GALBA + EvidenceModeler",
                        "github.com/ebp-nor/GenomeAnnotation"),
    "GCA_977110915.1": ("EBP-Nor pipeline: miniprot + GALBA + EvidenceModeler",
                        "github.com/ebp-nor/GenomeAnnotation"),
    "GCA_977091265.1": ("EBP-Nor pipeline: miniprot + GALBA + EvidenceModeler",
                        "github.com/ebp-nor/GenomeAnnotation"),
    "GCA_977110975.1": ("EBP-Nor pipeline: miniprot + GALBA + EvidenceModeler",
                        "github.com/ebp-nor/GenomeAnnotation"),
    "GCA_019671135.1": ("Funannotate v1.8.1, default parameters", "Aldossari & Ishii 2021, Microbiol Resour Announc 10:e00692-21"),
    "GCA_025094155.1": ("JGI Annotation Pipeline", "JGI MycoCosm portal Lichy1"),
    "GCF_041956525.1": ("JGI Annotation Pipeline", "JGI MycoCosm portal Rhipu1"),
    "GCA_003550305.1": ("JGI Annotation Pipeline", "Morin et al. 2019, New Phytol 222:1584-1598"),
    "JGI:Morwol1_v1.0": ("JGI Annotation Pipeline", "JGI MycoCosm portal Morwol1"),
    "GCA_040209145.1": ("JGI Annotation Pipeline", "JGI MycoCosm portal PhyblU21_2"),
}

SUB_SHORT = {"Mucoromycotina": "Mucoro.", "Mortierellomycotina": "Mortierello.", "Glomeromycotina": "Glomero."}
CELL_NAME = {
    "A": "A (B3 x R+)", "B": "B (B2 x R-)", "C": "C (B3 x R-)", "D": "D (B2 x R+)",
}

# ---------------------------------------------------------------- Methods 正文
METHODS_EN = [
 ("h", "Benchmark species selection for annotation comparison"),
 ("p", "To assess the reliability of the annotations delivered here, we selected a benchmark subset of 20 species "
       "from the 202-strain cohort and compared each strain's own annotation against the pre-existing third-party "
       "annotation of the same assembly at the protein level. Of the 202 strains, 119 were annotated with protein "
       "evidence only (BRAKER2) and 83 with protein plus RNA-seq evidence (BRAKER3). Selecting benchmarks without "
       "regard to the evidence basis of the reference annotation would therefore confound method performance with "
       "evidence availability: a protein-only annotation compared against a transcript-supported reference will "
       "systematically miss weakly expressed, fast-evolving or lineage-specific genes, and this deficit would be "
       "misread as lower accuracy. We therefore used a 2 × 2 factorial design in which one axis is the annotation "
       "mode applied on our side (BRAKER3 vs BRAKER2) and the other is the evidence basis of the reference "
       "annotation (transcript-supported, R+, vs protein-only, R−)."),
 ("h", "Reference annotation evidence tiers"),
 ("p", "Each reference annotation was assigned to one of four tiers by inspecting primary sources (publication "
       "methods sections, INSDC and MycoCosm records, and pipeline repositories). R1: the annotation was "
       "demonstrably produced using RNA-seq evidence. R2: the annotation was produced by the JGI Annotation "
       "Pipeline, for which transcript assemblies are a standard input; this tier is assigned at pipeline level "
       "and per-genome confirmation is pending. R3: the annotation used no transcript evidence, being based on "
       "protein homology only, ab initio prediction, or default parameter settings. R?: the evidence basis could "
       "not be established. R1 and R2 were grouped as R+ and R3 as R−; species assigned R? were excluded. "
       "A full per-species provenance record, including the software and version recorded by each submitter, is "
       "provided in Supplementary Table Sx."),
 ("h", "Selection rules and their order of application"),
 ("p", "Rules were applied in three ordered tiers, and within each tier they are exhaustive: a strain failing any "
       "Tier 0 gate was excluded without exception; Tier 1 quotas were then applied to the surviving set; Tier 2 "
       "ranking was used only to order species that compete for the same design cell."),
 ("li", "Tier 0 — hard gates. (G1) The strain is a member of the 202-strain cohort, matched uniquely by assembly "
        "accession. (G2) A third-party protein annotation exists for that assembly. (G3) Genome completeness in "
        "mucoromycota_odb12 genome mode is ≥ 98.5%, and assembly contiguity is either N50 ≥ 800 kb or "
        "chromosome-level. (G4) The reference annotation can be assigned to an evidence tier R1–R3. (G5) The "
        "reference protein set is retrievable."),
 ("li", "Tier 1 — design quotas. (D1) The 20 species are allocated 8 / 5 / 4 / 3 to cells A / B / C / D. "
        "(D2) Subphylum allocation is 12 Mucoromycotina, 5 Mortierellomycotina and 3 Glomeromycotina, "
        "proportional to the cohort (123 : 57 : 22) with deliberate enrichment of Glomeromycotina, the subphylum "
        "with the most fragmented assemblies; at least two glomeromycotan species are placed in the priority "
        "group. (D3) The priority group contains only evidence-matched cells, in the ratio A : B = 6 : 4."),
 ("li", "Tier 2 — within-cell ranking, applied in the following order: (i) higher genome completeness; "
        "(ii) greater assembly contiguity; (iii) higher reference evidence tier (R1 over R2 over R3); "
        "(iv) availability of a peer-reviewed reference-annotation publication that records software versions; "
        "(v) RefSeq reference-genome status; (vi) availability of RNA-seq runs already processed in this study; "
        "(vii) breadth of prior reuse of the strain in the literature."),
 ("h", "Application and composition of the benchmark set"),
 ("p", "Applying the Tier 0 gates to the 202-strain cohort retained 143 strains with a third-party annotation, "
       "78 of which reached ≥ 98.5% genome completeness and 48 of which also satisfied the contiguity "
       "requirement; the final 20 species were drawn from these 48. The priority group comprises the two "
       "evidence-matched cells (A, n = 6; B, n = 4), so that the transcript-supported and the protein-only halves "
       "of our delivery are each evaluated against references built on the same evidence class. The reserve group "
       "comprises four cell-C and three cell-D species, used only to estimate the effect of evidence type "
       "(Δ = C − A for the gain attributable to RNA-seq; Δ = D − B for the loss when RNA-seq is absent), "
       "together with three further evidence-matched species that reinforce cells A and B. The composition of the "
       "benchmark set is given in Table X."),
 ("h", "Limitations of the benchmark set"),
 ("p", "Three limitations apply. First, none of the reference annotations in this cohort is manually curated at "
       "the gene-model level: the RefSeq category 'reference genome' is a classification rather than evidence of "
       "curation, and all 12 GCF records in the cohort have paired GenBank records from which the RefSeq "
       "annotation is propagated. Comparisons therefore measure agreement between annotations and not accuracy "
       "against a truth standard. Second, four of the priority reference annotations were themselves produced by "
       "BRAKER-family pipelines (BRAKER 2.1.5, BRAKER1 v1.9, BRAKER v2.1.2, and BRAKER2 in protein-only mode), so "
       "agreement with them is partly a within-family comparison. Of the remaining six priority references, three use "
       "the JGI Annotation Pipeline and are methodologically independent, while three are EBP-Nor pipeline "
       "predictions based on protein homology alone. Third, the protein-only matched cell (B) could "
       "only be populated with Mucoromycotina species, because strains that combine a high-quality genome, a "
       "protein-only reference annotation and protein-only annotation on our side occur only in that subphylum; "
       "conclusions drawn from cell B are therefore restricted to Mucoromycotina."),
 ("h", "Table X. Benchmark species used for annotation comparison."),
 ("p", "Group, priority or reserve. Cell, design cell (A, BRAKER3 × R+; B, BRAKER2 × R−; C, BRAKER3 × R−; D, BRAKER2 × R+). Own mode, "
       "annotation mode applied in this study. Ref tier, evidence tier of the reference annotation (R1, RNA-seq "
       "evidence; R2, JGI Annotation Pipeline with transcript input; R3, no transcript evidence). N50, scaffold "
       "N50 of the assembly. BUSCO C, completeness in mucoromycota_odb12 genome mode. Reference proteins, "
       "number of protein sequences in the reference annotation. Reserve-group species are used only to estimate the effect of evidence type and are not "
       "included in the primary agreement statistics."),
]

RULES_MD = [
    ("Tier 0　硬性门槛（任一不满足即出局，无例外）", [
        "G1　该菌株属于本研究 202 株队列，以 Assembly accession 唯一锚定。",
        "G2　该装配存在第三方既有蛋白注释（表1「原始蛋白注释 = 有」）。",
        "G3　基因组质量：mucoromycota_odb12 基因组模式完整度 C ≥ 98.5%，且 N50 ≥ 800 kb 或装配为染色体级。",
        "G4　参考注释可归入证据等级 R1–R3（R? 一律出局）。",
        "G5　参考蛋白序列可取得（NCBI 可下载，或 JGI 门户可取）。",
    ]),
    ("Tier 1　设计配额（决定谁能进 20 株）", [
        "D1　设计格配额 A : B : C : D = 8 : 5 : 4 : 3。",
        "D2　亚门配额 Mucoromycotina 12 : Mortierellomycotina 5 : Glomeromycotina 3；"
        "Glomeromycotina 优先组不少于 2 株。",
        "D3　优先组只收证据匹配格，A : B = 6 : 4。",
    ]),
    ("Tier 2　同格内排序（仅用于决定同格内谁优先）", [
        "T1　基因组完整度（BUSCO C）高者优先。",
        "T2　装配连续性（N50）高者、或染色体级者优先。",
        "T3　参考注释证据等级高者优先（R1 > R2 > R3）。",
        "T4　参考注释有同行评审论文且记录了软件版本者优先。",
        "T5　RefSeq（GCF_）参考基因组优先。",
        "T6　队列内已有可复用 RNA-seq run 者优先。",
        "T7　菌株在文献中被广泛复用者优先。",
    ]),
]

BULLETS_CN = [
    "设计：2×2 因子——本方注释方式（BRAKER3 / BRAKER2）× 参考注释证据基础（R+ 含转录证据 / R− 仅蛋白）。"
    "202 株中 BRAKER2 119 株、BRAKER3 83 株；不作此区分会把「证据更薄」误读成「方法更差」。",
    "优先组 10 株全部来自两个证据匹配格（A 6 株、B 4 株），因此本方「有 RNA-seq 的一半」与"
    "「只有蛋白的一半」各自在公平条件下评估。",
    "备选组 10 株 = 错配格 C 4 株 + 错配格 D 3 株 + 同格补充 3 株；C、D 只用于量化证据类型效应"
    "（Δ = C − A 为多一套 RNA-seq 的增益，Δ = D − B 为少一套的损失），不并入主数字。",
    "筛选漏斗：202 → 143（有既有注释）→ 78（BUSCO C ≥ 98.5%）→ 48（且 N50 ≥ 800 kb 或染色体级）→ 20。",
    "亚门配额 12 : 5 : 3，对应队列 123 : 57 : 22；Glomeromycotina 有意加权，因其装配最破碎、最易被质疑。",
]


def build():
    from docx import Document
    from docx.shared import Pt, Cm, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.section import WD_ORIENT
    from docx.oxml.ns import qn

    d = pd.read_csv(CSV)
    doc = Document()

    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(9)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    for i in range(1, 5):
        h = doc.styles[f"Heading {i}"]
        h.font.name = "Times New Roman"
        h.element.rPr.rFonts.set(qn("w:eastAsia"), "黑体")
        h.font.color.rgb = RGBColor(0, 0, 0)
        h.font.size = Pt({1: 15, 2: 12.5, 3: 10.5, 4: 9.5}[i])
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = Cm(29.7), Cm(21.0)
    sec.left_margin = sec.right_margin = Cm(1.6)
    sec.top_margin = sec.bottom_margin = Cm(1.5)

    def para(t, bold=False, size=9, align=None, after=5, italic=False):
        p = doc.add_paragraph()
        r = p.add_run(t); r.bold = bold; r.italic = italic
        r.font.size = Pt(size)
        if align == "c":
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_after = Pt(after)
        return p

    def bullet(t, size=8.5):
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(t).font.size = Pt(size)
        p.paragraph_format.space_after = Pt(1)
        return p

    def table(headers, data, widths=None, font=7):
        t = doc.add_table(rows=1, cols=len(headers)); t.style = "Table Grid"
        for i, hh in enumerate(headers):
            c = t.rows[0].cells[i]; c.text = ""
            r = c.paragraphs[0].add_run(hh); r.bold = True; r.font.size = Pt(font)
            r.font.name = "Times New Roman"; r.element.rPr.rFonts.set(qn("w:eastAsia"), "黑体")
        for row in data:
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = ""
                r = cells[i].paragraphs[0].add_run(str(v)); r.font.size = Pt(font)
                r.font.name = "Times New Roman"; r.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
        if widths:
            for row in t.rows:
                for i, w in enumerate(widths):
                    row.cells[i].width = Cm(w)
        doc.add_paragraph().paragraph_format.space_after = Pt(1)
        return t

    # ---------- 标题 ----------
    doc.add_heading("基准物种遴选 · 最终简洁版", level=0)
    para("Mucoromycota 202 基因组注释项目 ｜ 批次一 · 任务 1.1 ｜ 归属：Technical Validation"
         "（规则部分亦可置于 Methods） ｜ 2026-09-14", size=9, align="c", after=10)

    doc.add_heading("一　要点", level=2)
    for b in BULLETS_CN:
        bullet(b)

    # ---------- 20 株表 ----------
    doc.add_heading("二　最终 20 株名单（优先 10 + 备选 10）", level=2)
    rows = []
    for _, r in d.iterrows():
        prov, ref = REF_PROV.get(r["Assembly_accession"], ("—", "—"))
        rows.append([
            r["组别"], r["标签"], CELL_NAME[r["设计格"]], r["物种"], r["菌株"],
            SUB_SHORT.get(r["亚门"], r["亚门"]), r["Assembly_accession"], r["本方注释方式"],
            r["参考注释证据等级"].split("（")[0], r["装配水平"], r["N50_kb"], r["基因组BUSCO_C_pct"],
            prov, r["既有注释基因数"],
        ])
    table(["组", "标", "设计格", "物种", "菌株", "亚门", "Assembly accession", "本方",
           "等级", "装配", "N50(kb)", "BUSCO C(%)", "参考注释流程", "基因数"],
          rows,
          widths=[0.8, 0.7, 1.7, 3.4, 1.8, 1.6, 2.6, 1.3, 0.8, 1.6, 1.1, 1.2, 4.3, 1.0], font=6.5)
    para("组：优先 / 备选。设计格 A = 本方 BRAKER3 × 参考 R+；B = 本方 BRAKER2 × 参考 R−；"
         "C = 本方 BRAKER3 × 参考 R−；D = 本方 BRAKER2 × 参考 R+。R+ = 含转录证据（R1 确证使用 RNA-seq，"
         "R2 JGI 流程以转录组为标准输入、逐株待核）；R− = 无转录证据。基因数为 NCBI 报告的基因总数"
         "（含 ncRNA/tRNA）；JGI 门户条目 NCBI 无记录，标 n.a.。", size=7.5)

    # ---------- 规则与优先级 ----------
    doc.add_heading("三　筛选规则与优先级（三级，按序应用）", level=2)
    for title, items in RULES_MD:
        doc.add_heading(title, level=3)
        for it in items:
            bullet(it, size=8.5)
    para("说明：Tier 0 是硬门槛，任一不满足即出局，无例外；只有通过 Tier 0 的物种才进入 Tier 1 配额分配；"
         "Tier 2 仅用于决定同一设计格内部的先后，不改变配额。", size=8)

    # ---------- 英文方法学原文 ----------
    doc.add_heading("四　可直接写入 Scientific Data 的英文原文", level=2)
    para("以下文字可整段粘贴。建议将「Benchmark species selection for annotation comparison」"
         "至「Selection rules and their order of application」置于 Methods 的注释设计小节，"
         "「Application and composition」与「Limitations」置于 Technical Validation 的开头。",
         size=8.5)
    for kind, text in METHODS_EN:
        if kind == "h":
            doc.add_heading(text, level=3)
        elif kind == "li":
            bullet(text, size=8.5)
        else:
            para(text, size=8.5, after=5)

    # ---------- 占位符 ----------
    doc.add_heading("五　投稿前需替换的占位符", level=2)
    table(["占位符", "需要填什么"],
          [["Table X", "本文正文中 20 株基准物种表的位置（若并入主表，请改用相应编号）"],
           ["Supplementary Table Sx", "逐株参考注释溯源表（流程、软件版本、证据类型、文献），建议随附表提交"],
           ["[repository]", "数据仓储 DOI 或访问号（本报告未使用，若正文需要请补）"]],
          widths=[4.5, 24.0], font=8)

    doc.save(DOCX)
    return DOCX


if __name__ == "__main__":
    print("written:", build())
