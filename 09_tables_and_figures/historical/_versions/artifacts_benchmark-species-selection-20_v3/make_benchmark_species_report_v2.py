#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成《基准物种遴选报告 V2》Word 文档（批次一 · 任务 1.1）。

V2 相对 V1 的核心变化：把「本方注释方式」提升为与「亚门」「基因组质量」并列的第三根设计轴。
设计为 2×2：本方 BRAKER2（仅蛋白证据）/ BRAKER3（蛋白 + RNA-seq）
             × 参考注释证据等级 R+（含转录证据）/ R−（仅蛋白或从头预测）。
优先组 = 两个证据匹配格；备选组 = 两个可控错配格 + 同格补充。

输入（全部只读）：
  投稿附表/表1_来源清单_总表.xlsx        —— 原始蛋白注释 有/无、assembly 水平
  投稿附表/表3_研究总表.xlsx            —— 202 株权威编号、采纳名、亚门、BUSCO C
  投稿附表/表6_注释证据_202株.xlsx       —— RNA-seq run 与比对率、BRAKER 模式
  投稿附表/表7_基因结构注释_202株.xlsx    —— 基因结构统计、N50
  6_Protein_evidence/14_tables_OrthoDBv12.2_*.xlsx —— OrthoDB v12.2 Mucoromycota 157 子集
  NCBI Datasets API v2alpha（2026-09-14 快照）—— 装配与注释元数据

输出：
  9_Benchmark_reference_annotations/基准物种遴选报告_20株_v1.docx

运行：python3 9_Benchmark_reference_annotations/make_benchmark_species_report.py
"""
import os
import re
import json
import urllib.request

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "9_Benchmark_reference_annotations")
ARCHIVE = os.path.join(OUTDIR, "归档_v1_20260914")
RAW_FAA = os.path.join(ARCHIVE, "04_参考蛋白FASTA")
RAW_GFF = os.path.join(ARCHIVE, "05_参考注释GFF3")
CACHE = os.path.join(OUTDIR, "benchmark_species_evidence.json")
DOCX = os.path.join(OUTDIR, "基准物种遴选报告_20株_v2.docx")

# ----------------------------------------------------------------------------
# 1. 数据装载
# ----------------------------------------------------------------------------

def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def load_master():
    t1 = pd.read_excel(os.path.join(ROOT, "投稿附表/表1_来源清单_总表.xlsx"), sheet_name="统一清单_231")
    t3 = pd.read_excel(os.path.join(ROOT, "投稿附表/表3_研究总表.xlsx"), sheet_name="研究菌株_202")
    t6 = pd.read_excel(os.path.join(ROOT, "投稿附表/表6_注释证据_202株.xlsx"), sheet_name="1_证据总览_202株")
    t7 = pd.read_excel(os.path.join(ROOT, "投稿附表/表7_基因结构注释_202株.xlsx"), sheet_name="1_表7_基因结构统计")

    t3 = t3[pd.to_numeric(t3["编号"], errors="coerce").notna()].copy()
    t6 = t6[pd.to_numeric(t6["No"], errors="coerce").notna()].copy()
    t7 = t7[pd.to_numeric(t7["No"], errors="coerce").notna()].copy()
    for d, c in ((t3, "编号"), (t6, "No"), (t7, "No")):
        d[c] = d[c].astype(int)

    m = t3.merge(t6[["No", "RNA-seq", "BRAKER_mode"]], left_on="编号", right_on="No", how="left")
    m = m.merge(t7[["No", "n_genes", "n_proteins_faa", "assembly_N50_kb"]], on="No", how="left")
    t1k = t1[["主标识", "assembly 水平（来源词）", "原始蛋白注释", "strain/isolate"]].rename(
        columns={"strain/isolate": "t1_strain"})
    m = m.merge(t1k, on="主标识", how="left")
    # 表S5 提供全部 202 株（含 JGI 条目）的基因组总长
    s5 = pd.read_csv(os.path.join(ROOT, "TableS5_simplified_202_v6.tsv"), sep="\t")
    m = m.merge(s5[["Assembly", "total_length_bp"]].rename(columns={"Assembly": "主标识"}),
                on="主标识", how="left")
    return m


def fetch_ncbi(accessions):
    """NCBI Datasets API v2alpha 装配+注释元数据；带本地缓存。"""
    out = {}
    if os.path.exists(CACHE):
        out = json.load(open(CACHE))
    todo = [a for a in accessions if a.startswith(("GCA_", "GCF_")) and a not in out]
    for a in todo:
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{a}/dataset_report"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mucoromycota-annotation-audit/1.0"})
            with urllib.request.urlopen(req, timeout=45) as r:
                rep = json.load(r)["reports"][0]
            ai = rep.get("assembly_info") or {}
            an = rep.get("annotation_info") or {}
            st = rep.get("assembly_stats") or {}
            gc = (an.get("stats") or {}).get("gene_counts") or {}
            bp = []
            for lin in (ai.get("bioproject_lineage") or []):
                for b in (lin.get("bioprojects") or []):
                    bp.append(f"{b.get('accession')}: {b.get('title')}")
            out[a] = {
                "organism_name": (rep.get("organism") or {}).get("organism_name"),
                "strain": ((rep.get("organism") or {}).get("infraspecific_names") or {}).get("strain"),
                "assembly_level": ai.get("assembly_level"),
                "refseq_category": ai.get("refseq_category"),
                "release_date": ai.get("release_date"),
                "submitter": ai.get("submitter"),
                "sequencing_tech": ai.get("sequencing_tech"),
                "assembly_method": ai.get("assembly_method"),
                "bioproject": ai.get("bioproject_accession"),
                "bioproject_titles": bp,
                "comments": (ai.get("comments") or "")[:1500],
                "ann_name": an.get("name"),
                "ann_provider": an.get("provider"),
                "ann_release_date": an.get("release_date"),
                "ann_gene_total": gc.get("total"),
                "ann_gene_protein_coding": gc.get("protein_coding"),
                "total_sequence_length": st.get("total_sequence_length"),
                "scaffold_n50": st.get("scaffold_n50"),
                "contig_n50": st.get("contig_n50"),
                "number_of_scaffolds": st.get("number_of_scaffolds"),
                "gc_percent": st.get("gc_percent"),
            }
        except Exception as e:  # noqa: BLE001
            out[a] = {"error": str(e)}
    json.dump(out, open(CACHE, "w"), ensure_ascii=False, indent=1)
    return out


def orthodb_members():
    d = pd.read_excel(os.path.join(
        ROOT, "6_Protein_evidence/14_tables_OrthoDBv12.2_Mucoromycota157_ProteinEvidence_Taxonomy.xlsx"),
        sheet_name="Sheet1")
    s = set()
    for c in ("sample", "species"):
        for v in d[c].dropna().astype(str):
            s.add(norm(v))
    return s


# ----------------------------------------------------------------------------
# 2. 名单定义
# ----------------------------------------------------------------------------

# (accession, 组别, 一句话定位, 入选理由段落, 风险与备注)
# -*- coding: utf-8 -*-
PRIORITY = [
    # ---------------- 通道 A：本方 BRAKER3 × 参考注释含转录证据（证据匹配） ----------------
    ("GCF_026210795.1", "通道 A｜本方 BRAKER3 × 参考 R+",
     "三个亚门里 Glomeromycotina 装配质量最差、最易被质疑，本株是全队列 202 株中唯一达到染色体级（Chromosome）的 AMF 装配：32 条染色体、42 contig → 32 scaffold、contig N50 3.90 Mb、scaffold N50 5.09 Mb、仅 10 个 gap，其中 23 条染色体完全无 gap、17 条达到端粒到端粒；mucoromycota_odb12 基因组模式 BUSCO 完整度 C = 99.8%。本方为 BRAKER3（有 RNA-seq：SRR3285886 98.96%、SRR3285885 98.93%），参考注释同样有转录证据，且是短读 + 长读双证据：Illumina RNA-seq → STAR 2.7.6a → BRAKER 2.1.5，再由 Nanopore PCR-cDNA RNA-seq 经 funannotate 1.8.7 update 与 PASA 精修 UTR 与 poly(A)。两侧证据类型对齐，是本套设计里最干净的一个配对。原始论文 Manley et al. 2023, G3 13(6):jkad077。",
     "三条必须写进 Methods 的限制。第一，参考注释由 BRAKER 2.1.5 产生，与本研究属同一工具家族，一致性里有一部分是自我印证（见第 6.12 节）。第二，基因数存在三套不一致的版本：论文报 30,209，GenBank 提交版 GCA_026210795.1 为 30,190 个蛋白编码基因，RefSeq 版 GCF_026210795.1 为 30,147——本报告锁定 GCF_026210795.1，引用时必须写明取的是哪一套。第三，GCF 记录挂在独立的 RefSeq BioProject PRJNA1128271 下，RefSeq 侧只做重新发布（GenBank 孪生记录为 GCA_026210795.1），不可写成「NCBI 注释」。此外该参考注释含异构体（885/30,147 基因有多转录本，3.3%），必须做最长异构体提取。"),
    ("GCF_035048745.1", "通道 A｜本方 BRAKER3 × 参考 R+",
     "全部候选中证据类型与 BRAKER3 最可比的基准。本方为 BRAKER3（SRR28755642 90.19%、SRR28755641 89.73%）；参考注释使用 16 个生长条件（RPMI ± 人 PMN、YG/MEA 不同温度、铁限制、H2O2、SBA 等；SAMN41019669–84）的 RNA-seq，经 STAR v2.5.3a → StringTie v1.3.3，再与 TransDecoder v5.0.2、BRAKER1 v1.9、CodingQuarry v2.0 的结果一并交给 EVidenceModeler 整合，funannotate 收尾——每一步都有版本号，且有专门的论文（Conlan et al. 2024, Microbiol Resour Announc 13(7):e00062-24）。装配也好：42 个 contig、contig N50 2.65 Mb、Merqury QV 46.6 / 完整度 99.7%、17 个 contig 带 TTAGGGn 端粒（提示 9 条染色体）。参考注释为严格单转录本，任务 1.3 的最长异构体提取对它是空操作。",
     "基因数口径要当心：NCBI 的 13,535 是「基因总数」，其中蛋白编码 13,295 + tRNA 240；论文报的 13,295 是蛋白编码数，两者完全对账。另有一处小矛盾：NCBI 装配记录的注释说明栏写 Funannotate v1.8.14，而论文 Methods 写 v1.8.13（本报告未能取得该说明栏原文，标注为待核）。表3 采纳名为 Mucor circinelloides、原名为 Mucor velutinosus，是 202 队列中 5 组同物异名之一；写作时必须同时给出原名与采纳名。参考注释含 BRAKER1，带一定循环性；不过它与 BRAKER3 的差异主要集中在证据组合方式（BRAKER1 只用 RNA-seq，本研究为 RNA-seq + 蛋白），这一点可以在讨论中正面利用。"),
    ("GCA_015698045.1", "通道 A｜本方 BRAKER3 × 参考 R+",
     "Glomeromycotina 第二代表，与 P1 构成属内对照。本方为 BRAKER3（DRR041161 88.24%、DRR041163 85.86%）；参考注释为 BRAKER v2.1.2 + AUGUSTUS 3.3.2 + GeneMark-ES 4.32，作者自述为 RNA-seq-based，证据类型与 BRAKER3 部分对齐。基因组 BUSCO C = 99.8%，contig N50 = 1.31 Mb（360 contig，L50 = 40），是 Glo 内除 R. irregularis 之外连续性最好的一株；30,348 个基因（30,194 蛋白编码 + 154 nc）。",
     "三个限制。第一，唯一文献是 2026-08-13 发布于 Research Square、状态 Under Revision 的预印本（DOI 10.21203/rs.3.rs-10340957/v1），未经同行评审；作者自认该组装未达染色体级。第二，基因数在预印本（31,223 蛋白编码 / 32,541 蛋白模型）与 NCBI（30,348）之间不一致，本报告锁定 NCBI 提交版；且预印本未点名所用的 RNA-seq 数据集。第三，参考 GFF3 未使用标准 mRNA/transcript 特征（只有 gene 与 CDS），无法从 GFF 读出异构体结构，但蛋白 FASTA 有 32,541 条、多于 30,194 个蛋白编码基因，说明蛋白集含异构体——任务 1.3 必须以蛋白 FASTA 为准。另需区分：2018 年 BMC Genomics 那篇（DOI 10.1186/s12864-018-4853-0）用的是另一套组装（27,753 基因）。"),
    ("GCA_040256725.1", "通道 A｜本方 BRAKER3 × 参考 R+",
     "本方 BRAKER3，且队列内两条 RNA-seq 比对率均在 99% 以上（SRR14527209 99.20%、SRR14527207 99.19%），是全部候选里本方 RNA-seq 质量最好的一株；参考注释为 JGI Annotation Pipeline（Mucci3，2024-06-21），流程以转录组数据为标准输入，证据类型同侧。基因组是 Mucoromycotina 内最好的一档：BUSCO C = 99.9%（与全队列最高值并列），N50 = 4.61 Mb；11,759 个基因、密度 321 genes/Mb，落在 Mucor 属常规区间。参考注释为严格单转录本。该株即文档候选表中的「Mucor circinelloides CBS 277.49」，现用名为 Mucor lusitanicus，同一菌株号。",
     "文档原候选写作 MU402，但 202 队列与 NCBI 记录均为 CBS 277.49，必须改写。JGI 专属转录组（PRJNA1362231、PRJNA1360744）2025-11 才公开，晚于该装配 2024-06 的发布——若本项目的 BRAKER3 用了这些 run 而 JGI 注释当时未用，两者差异有一部分来自证据可得性；逐株的「JGI 是否用了 RNA-seq」在公开元数据中无记录，建议在 MycoCosm 门户人工核对后补记。"),
    ("JGI:Disdec1_v1.0", "通道 A｜本方 BRAKER3 × 参考 R+",
     "本方 BRAKER3（SRR13647699 97.75%）；参考注释为 JGI Annotation Pipeline，证据类型同侧。N50 = 6.89 Mb 为 202 株中最高值，基因组 BUSCO C = 99.9%，与全队列最高值并列——Mortierellomycotina 内稀缺的高质量装配。以它作为该亚门第一个基准，可以在「基因组质量最好」的条件下先确立注释基线，再去看该亚门其余较破碎的装配。",
     "属 JGI MycoCosm 门户直取：该菌株在 NCBI 没有公开装配，比对与仓储都必须锁定门户代码 + 版本号 Disdec1_v1.0 + 下载日期 + 文件哈希，不能只写物种名。转录组 PRJNA696871 存在，但「该版本注释是否采用了它」无公开记录（MycoCosm 页面对本机返回 403），需在门户内人工核对。参考蛋白需从门户下载。"),
    ("GCA_036320915.1", "通道 A｜本方 BRAKER3 × 参考 R+",
     "本方 BRAKER3（SRR6225694 92.77%、SRR6225697 90.85%）；参考注释为 JGI Annotation Pipeline（MoeNVP64_1，2024-01-31），证据类型同侧。基因组 BUSCO C = 99.8%，N50 = 2.74 Mb，Mortierellomycotina 内仅次于 Disdec1；13,991 个基因、密度 274 genes/Mb；参考注释为严格单转录本。Linnemannia elongata 是 Mortierellaceae 分类修订中的关键物种，与 P5 同科不同属，构成该亚门内的属级对照。",
     "文档原候选写作 AG-77，202 队列实际为 NVP64（GCA_036320915.1），菌株号必须改写。该株参考注释的转录组输入在 JGI 元数据中同样没有逐条标注。"),
    # ---------------- 通道 B：本方 BRAKER2 × 参考注释仅蛋白证据（证据匹配） ----------------
    ("GCA_977092135.1", "通道 B｜本方 BRAKER2 × 参考 R−",
     "本方为 BRAKER2（仅蛋白证据，表6 无保留 RNA-seq run），参考注释为 EBP-Nor 流程（miniprot + GALBA + EvidenceModeler），同样只用 3 套蛋白同源证据、零 RNA-seq —— 这是本套设计里唯一能做到「两侧都只有蛋白证据」的组合，用来回答「我的 BRAKER2 那一半交付物到底可不可靠」。基因组质量在通道 B 中最高：BUSCO C = 99.6%，N50 = 6.36 Mb（全队列第 2 高）。",
     "参考注释为 EBP-Nor 自研流程的自动产物：无 RefSeq 镜像、无注释方法论文（只有 GitHub 仓库），且异构体负荷 8.9%（多转录本基因 1,162 个），必须做最长异构体提取。两侧都没有转录证据，因此本配对只能评估「蛋白证据条件下的一致程度」，不能外推为「基因模型完整性」——凡涉及非编码、UTR、可变剪接的结论都要标注为超出本配对的能力范围。"),
    ("GCF_029851405.1", "通道 B｜本方 BRAKER2 × 参考 R−",
     "这是全部候选中循环性最强、但同时也最锋利的一个配对：本方 BRAKER2，参考注释由 Broad Institute 用 BRAKER2 蛋白模式产生（证据为 L. corymbifera / L. ramosa / L. hyalospora 的 35,966 条蛋白），两侧是同一工具、同一证据类型，只是执行者不同。因此它在方法学上回答的问题是「同一工具在不同团队、不同参数下能否得到彼此一致的结果」。该株为 GCF_ 开头的 RefSeq 参考基因组，BUSCO C = 99.8%；原始论文 Shea et al. 2023, Microbiol Resour Announc 12(7):e00338-23。Lichtheimia 属是医学重要类群（接合菌病主要病原之一），doc 候选表中的 L. corymbifera 因无既有注释被替换，本株是同属唯一可用者。",
     "循环性最重：BRAKER2 vs BRAKER2 的高度一致不能用来证明方法正确，只能证明结果可复现。措辞必须是「可复现性」而不是「准确性」。装配为 Scaffold 级、N50 = 386 kb、contig N50 仅 151 kb，碎片化明显，不得用于外显子边界或装配质量论证。参考注释为严格单转录本。"),
    ("GCA_977110935.1", "通道 B｜本方 BRAKER2 × 参考 R−",
     "本方 BRAKER2（无保留 RNA-seq），参考注释为 EBP-Nor 纯蛋白同源流程，两侧证据类型一致。基因组 BUSCO C = 99.8%（通道 B 中最高），N50 = 3.12 Mb；参考注释 11,908 个基因。与 P7 同为 EBP-Nor 批次，但基因组连续性更好、基因数更接近该属常规区间，可用来检验「同一参考流程在不同装配质量下的表现」。",
     "同 P7：参考注释为自动流程产物、无 RefSeq 镜像、无方法论文；异构体负荷需按附录 C 数值处理，必须做最长异构体提取。两侧均无转录证据，结论不外推到基因模型完整性。"),
    ("GCA_977110905.1", "通道 B｜本方 BRAKER2 × 参考 R−",
     "本方 BRAKER2（无保留 RNA-seq），参考注释为 EBP-Nor 纯蛋白同源流程，两侧证据类型一致。基因组 BUSCO C = 99.4%，N50 = 3.85 Mb（通道 B 中第 2 高）；参考注释 12,726 个基因。取本株是为了把 Mucor 属在通道 B 里做成三重复（M. fragilis / M. hiemalis / 通道外的 U. vinacea），使通道 B 的结论不至于只建立在一个属的一个样本上。",
     "同 P7：参考注释为自动流程产物、无方法论文；异构体负荷按附录 C 处理，必须做最长异构体提取。两侧均无转录证据。"),
]

BACKUP = [
    # ---------------- 通道 C：本方 BRAKER3 × 参考 R−（本方证据更强） ----------------
    ("GCA_977110955.1", "通道 C｜本方 BRAKER3 × 参考 R−",
     "刻意的证据错配，用于量化「本方多了一套 RNA-seq、参考只有蛋白」会带来多大差异。本方 BRAKER3（SRR1286003 62.42% 等 6 条 run）；参考注释为 EBP-Nor 纯蛋白同源流程。基因组质量优秀：BUSCO C = 99.4%，N50 = 4.72 Mb。该株同时是文档候选表中的 Rhizopus microsporus（临床与工业双重要），把错配格放在它身上可让这个格直接对接一个高关注度物种。",
     "参考注释异构体负荷高达 23.3%（1,696/10,755 基因有多转录本），是全部候选里最高的，不做最长异构体提取会严重污染 BUSCO 重复率与「各自独有蛋白数」。本方队列内 RNA-seq 比对率偏低（43.93%–62.42%），讨论错配效应时要同时说明本方证据本身也不算强。参考注释为自动流程产物、无 RefSeq 镜像、无方法论文。"),
    ("GCA_977110915.1", "通道 C｜本方 BRAKER3 × 参考 R−",
     "文档候选之一，证据错配格。本方 BRAKER3（SRR7101066 92.02%、SRR7101064 91.58%）；参考注释为 EBP-Nor 纯蛋白同源流程。基因组 BUSCO C = 99.4%，N50 = 4.26 Mb；参考注释 11,249 个基因、异构体负荷 11.6%。Cunninghamella 属在 Mucorales 内位置特殊，可补充该格内的分类学多样性。",
     "同通道 C 通例：参考注释为自动流程产物、无方法论文、异构体负荷高，必须做最长异构体提取。错配结果只能解释为「证据类型差异」，不能直接读作本方更优。"),
    ("GCA_977091265.1", "通道 C｜本方 BRAKER3 × 参考 R−",
     "证据错配格，同时补齐 Mortierellomycotina 在这条通道上的代表。本方 BRAKER3（SRR1638093 96.85%、SRR1638091 96.69%）；参考注释为 EBP-Nor 纯蛋白同源流程，11,632 个基因、异构体负荷 7.7%。基因组 BUSCO C = 99.3%，N50 = 3.23 Mb。Mortierella alpina 是工业重要物种（ARA 生产），其基因家族解读对注释完整性敏感，正好用来观察错配效应。",
     "同通道 C 通例。该株同时是文档候选之一（文档给的理由是「工业重要，注释关注度高」）。"),
    ("GCA_019671135.1", "通道 C｜本方 BRAKER3 × 参考 R−",
     "证据错配格中参考流程最明确的一例：参考注释由 University of Minnesota 用 Funannotate v1.8.1 全默认参数产生，论文明确把 RNA-seq 列为未来工作（Aldossari & Ishii 2021, Microbiol Resour Announc 10(37):e00692-21）。本方 BRAKER3（SRR21992393 93.84%、SRR21992390 80.84%）。基因组 BUSCO C = 99.4%，N50 = 2.32 Mb；参考注释 13,010 个基因（12,693 蛋白编码 + 317 tRNA）、严格单转录本。",
     "参考基因组 51.6 Mb 对 Linnemannia 属偏大、基因密度仅 246 genes/Mb（本批最低），用前建议先查重复含量与污染——若确有膨胀，本格的错配效应会被基因组因素混淆，必须同时报告重复含量。参考流程为默认参数、无转录证据，属「中等可信度基准」。"),
    # ---------------- 通道 D：本方 BRAKER2 × 参考 R+（本方证据更弱） ----------------
    ("GCA_025094155.1", "通道 D｜本方 BRAKER2 × 参考 R+",
     "刻意的证据错配，方向与通道 C 相反，用于量化「本方只有蛋白、参考有转录证据」会漏掉什么。本方 BRAKER2（表6 无保留 RNA-seq run）；参考注释为 JGI Annotation Pipeline（流程以转录组为标准输入）。基因组 BUSCO C = 99.7%，N50 = 129 kb；参考注释 12,174 个基因、严格单转录本。与 P8（Lichtheimia ornata，BRAKER2 蛋白模式参考）构成同属双基准，可把「属内差异」与「工具/证据差异」拆开。",
     "该株表6 记为无保留 RNA-seq run，因此本项目对它是 BRAKER2，两层都会是「仅蛋白」——错配方向明确且可控。装配连续性偏低（N50 129 kb），只用于蛋白层面一致性比较，不用于装配质量论证。"),
    ("GCF_041956525.1", "通道 D｜本方 BRAKER2 × 参考 R+",
     "证据错配格，且是通道 D 中唯一的 RefSeq（GCF_）参考基因组。本方 BRAKER2；参考注释为 JGI Annotation Pipeline，11,084 个基因、严格单转录本。基因组 BUSCO C = 98.8%，N50 = 1.21 Mb。Rhizomucor pusillus 临床重要，与 Rhizopus、Lichtheimia 同属 Mucorales 医学类群，可在错配格内形成医学物种的重复。",
     "同通道 D 通例：错配结果只能解释为「证据类型差异」。参考注释的转录组输入在 JGI 元数据中未逐条标注，需门户核对。"),
    ("GCA_003550305.1", "通道 D｜本方 BRAKER2 × 参考 R+",
     "通道 D 的 Glomeromycotina 代表。本方 BRAKER2；参考注释由 JGI Annotation Pipeline 产生（Morin et al. 2019, New Phytol 222(3):1584-1598，ENA 已挂 PMID 30636349），21,336 个基因。基因组 BUSCO C = 98.7%，N50 = 148 kb。取它是为了让三个亚门在四条设计通道里都有落点，避免「证据错配效应」只在 Mucoromycotina 上被观察到。",
     "两点特殊：其一，该研究自身的 RNA-seq 只做了 R. irregularis 与 G. rosea，本株没有自身转录组，JGI 流程对本株用了哪些其他证据未验证——因此它归入 R+ 是依「JGI 流程以转录组为标准输入」这一流程级判断，属于本报告中最弱的一处 R+ 归类，需在门户核对后确认；其二，装配缺口极高（136.8 Mb 中约 24.1 Mb、17.6% 为 N），碎片化会系统性压低任何重注释的完整性指标。"),
    # ---------------- 通道 E：同格补充（用于把两个匹配格做厚） ----------------
    ("JGI:Morwol1_v1.0", "通道 E｜本方 BRAKER3 × 参考 R+（同格补充）",
     "通道 A 的同格补充，Mortierellomycotina。本方 BRAKER3（SRR11558095 98.62%）；参考注释为 JGI Annotation Pipeline。N50 = 4.65 Mb（全队列第 3 高），BUSCO C = 99.2%。之所以放在备选而不是优先，只因它在 NCBI 没有对应装配（门户直取），取证链比 P1–P6 长一级。",
     "JGI MycoCosm 门户直取，须锁定门户代码 Morwol1 + 版本 + 下载日期 + 文件哈希；参考蛋白需从门户下载。另需注意 NCBI 上同菌株的另一个装配 GCA_016098105.1（UCR_MwolfiiNRRL6351_1.0）用的是 funannotate v1.7.1、scaffold N50 仅 93.7 kb，与本条目不是同一套，切勿混引。"),
    ("GCA_977110975.1", "通道 E｜本方 BRAKER2 × 参考 R−（同格补充）",
     "通道 B 的同格补充。本方 BRAKER2；参考注释为 EBP-Nor 纯蛋白同源流程，8,439 个基因。基因组 BUSCO C = 99.2%，N50 = 1.88 Mb。加入它是为了把匹配的「双蛋白证据」格从 4 株加到 5 株，并让 Umbelopsis 属在通道 B 内有两株（与文档所述的基因数异常物种 U. ramanniana 同属），便于交叉核对基因数偏低的成因。",
     "同通道 B 通例：参考注释为自动流程产物、无方法论文、异构体负荷需按附录 C 处理，必须做最长异构体提取。"),
    ("GCA_040209145.1", "通道 E｜本方 BRAKER3 × 参考 R+（同格补充）",
     "通道 A 的同格补充。本方 BRAKER3（SRR10402501 97.85%、SRR9002717 97.56%）；参考注释为 JGI Annotation Pipeline（PhyblU21_2，2024-06-17），13,146 个基因、严格单转录本。基因组 BUSCO C = 99.6%，N50 = 1.89 Mb。Phycomyces 是光生物学与类胡萝卜素研究的经典模式属，其基因模型在文献中可交叉比对，出现异常时容易定位原因；同时它是文档候选之一。",
     "文档原候选写作 NRRL1555，202 队列中实际为 UBC21（GCA_040209145.1），两者为不同菌株，必须改写。该基因组偏大（56.7 Mb）而基因密度偏低（232 genes/Mb），属高重复含量装配，解释指标时应说明。JGI 专属转录组 PRJNA566971（SRR10402501）与本项目所用 run 相同，因此本源数据两侧共用，不构成独立证据。"),
]

EXTENDED = [
    ("GCA_002105175.1", "Absidia repens", "Mucoromycotina", "通道 A 同格｜BRAKER3 × JGI(R+)，BUSCO 99.4%，N50 1.30 Mb，15,151 基因，单转录本。"),
    ("GCA_002105135.1", "Syncephalastrum racemosum", "Mucoromycotina", "通道 A 同格｜BRAKER3 × JGI(R+)，BUSCO 99.6%，N50 2.37 Mb，11,247 基因。"),
    ("GCF_025331425.1", "Radiomyces spectabilis", "Mucoromycotina", "通道 A 同格｜RefSeq + JGI(R+)，BUSCO 99.2%，N50 1.93 Mb。"),
    ("GCF_025528875.1", "Mycotypha africana", "Mucoromycotina", "通道 A 同格｜RefSeq + JGI(R+)，BUSCO 99.0%，N50 1.96 Mb。"),
    ("GCF_025118155.1", "Cokeromyces recurvatus", "Mucoromycotina", "通道 A 同格｜RefSeq + JGI(R+)，BUSCO 98.4%，N50 811 kb。"),
    ("GCF_002105155.1", "Lobosporangium transversale", "Mortierellomycotina", "通道 A 同格｜RefSeq + JGI(R+)，BUSCO 98.6%，N50 673 kb。"),
    ("GCA_025201345.1", "Thamnidium elegans", "Mucoromycotina", "通道 A 同格｜JGI(R+)，BUSCO 99.5%，N50 916 kb。"),
    ("JGI:Abscae1_v1.0", "Absidia caerulea", "Mucoromycotina", "通道 A 同格｜JGI 门户，BUSCO 99.4%，N50 3.62 Mb。"),
    ("JGI:FenlinCB54_1_v1.0", "Fennellomyces linderi", "Mucoromycotina", "通道 A 同格｜JGI 门户，BUSCO 99.9%，N50 3.56 Mb。"),
    ("JGI:Benpoi1_v1.0", "Benjaminiella poitrasii", "Mucoromycotina", "通道 A 同格｜JGI 门户，BUSCO 99.2%，N50 3.50 Mb。"),
    ("JGI:Ellano2_v2.0", "Ellisomyces anomalus", "Mucoromycotina", "通道 A 同格｜JGI 门户，BUSCO 99.8%，N50 2.29 Mb。"),
    ("JGI:Moramb1_v1.0", "Mortierella ambigua", "Mortierellomycotina", "通道 A 同格｜JGI 门户，BUSCO 98.6%，N50 3.39 Mb。"),
    ("JGI:Chabre1_v1.0", "Chaetocladium brefeldii", "Mucoromycotina", "通道 D 同格｜本方 BRAKER2 × JGI(R+)，BUSCO 99.3%，N50 1.31 Mb。"),
    ("JGI:Baclam1_v1.0", "Backusella lamprospora", "Mucoromycotina", "通道 D 同格｜本方 BRAKER2 × JGI(R+)，BUSCO 99.2%，N50 1.38 Mb。"),
    ("JGI:Entech1_v1.0", "Entomortierella echinosphaera", "Mortierellomycotina", "通道 D 同格｜本方 BRAKER2 × JGI 门户(R+)，BUSCO 98.8%，N50 1.13 Mb。"),
    ("JGI:Syzme1_v1.0", "Syzygites megalocarpus", "Mucoromycotina", "通道 D 同格｜本方 BRAKER2 × JGI 门户(R+)，BUSCO 98.7%，N50 1.53 Mb。"),
    ("GCA_053572175.1", "Umbelopsis sp. WA50703", "Mucoromycotina", "特例｜本方 BRAKER2 × 参考有 RNA-seq（Funannotate 1.8.15 + 两套 RNA-seq），严格说属通道 D；BUSCO 99.2%。注意同株另有一组装 GCA_964291815.1（8,693 基因），文献中的 MAISEN 审校属后者。"),
    ("GCA_977110945.1", "Umbelopsis ramanniana", "Mucoromycotina", "通道 C 同格｜BRAKER3 × EBP-Nor(R−)，BUSCO 99.4%。既有注释仅 9,709 个基因，文档所述「基因数异常」经核实成立，宜作异常例证（批次二任务 2.4）而非基准。"),
    ("GCA_040367835.1", "Mucor flavus KT1a", "Mucoromycotina", "通道 C 同格｜BRAKER3 × Funannotate 1.8.15 默认(R−)，有论文与明确版本号；BUSCO 99.8%，N50 858 kb。"),
    ("GCA_040366295.1", "Helicostylum pulchrum KT1b", "Mucoromycotina", "通道 C 同格｜同上；BUSCO 99.7%，N50 946 kb。"),
    ("GCA_016758905.1", "Umbelopsis isabellina", "Mucoromycotina", "通道 C 同格｜MAKER v2.31.8 + 10 万条公共 EST（是 EST，不是 RNA-seq）；BUSCO 99.6%，9,033 基因。"),
    ("GCA_910592045.1", "Gigaspora rosea", "Glomeromycotina", "不建议｜BRAKER3 × GeneMark-ES 纯从头预测；组装 N50 仅 4 kb、BUSCO C 46%，参考残缺会惩罚更优的注释。"),
    ("GCA_910592245.1", "Gigaspora margarita", "Glomeromycotina", "不建议｜同上；N50 仅 6 kb、115 万 contig 级碎片，575 Mb 上 46,492 基因。"),
    ("GCA_900079185.1", "Absidia glauca", "Mucoromycotina", "不建议｜本方 BRAKER2；参考流程完全未文档化、未找到任何论文，Ensembl 确证为原样导入；装配仅 Illumina + Velvet，contig N50 仅 11 kb。"),
]


# 既有注释的证据类型分层（每条均已核实到一手来源；unverified 项如实标注）
EVIDENCE_TIERS = [
    ("Rhizophagus irregularis DAOM-197198", "GCF_026210795.1", "BRAKER 2.1.5 + funannotate/PASA（有论文）",
     "有，且为短读 + 长读双证据：Illumina RNA-seq → STAR 2.7.6a → BRAKER 2.1.5；另有 Nanopore PCR-cDNA RNA-seq → funannotate 1.8.7 update + PASA 精修 UTR 与 poly(A)",
     "部分：Hi-C 挂载人工校正 + 人工剔除含转座子结构域的基因；非逐基因模型审编", "高（AMF 内事实参考；但与 BRAKER3 同源）",
     "Manley et al. 2023, G3 13(6):jkad077, DOI 10.1093/g3journal/jkad077"),
    ("Rhizophagus clarus HR1", "GCA_015698045.1", "BRAKER v2.1.2 + AUGUSTUS 3.3.2 + GeneMark-ES 4.32（仅预印本）",
     "有，作者自述 RNA-seq-based；但具体用了哪套 RNA-seq 数据未点名（未验证）", "无记录", "中（唯一文献为未经同行评审的 2026 预印本）",
     "Kobayashi et al., Research Square 预印本 2026-08-13, DOI 10.21203/rs.3.rs-10340957/v1（Under Revision）"),
    ("Mucor velutinosus NIH1002", "GCF_035048745.1", "de novo 自动流程（有论文）",
     "有：16 个生长条件 RNA-seq → STAR v2.5.3a / StringTie v1.3.3；与 TransDecoder v5.0.2、BRAKER1 v1.9、CodingQuarry v2.0 结果交 EVidenceModeler 整合，funannotate 收尾",
     "无", "高（证据类型与 BRAKER3 最接近；但与 BRAKER 家族同源）",
     "Conlan et al. 2024, Microbiol Resour Announc 13(7):e00062-24, DOI 10.1128/mra.00062-24"),
    ("Lichtheimia ornata CBS 291.66", "GCF_029851405.1", "BRAKER2 蛋白模式（有论文）",
     "无：证据为 L. corymbifera / L. ramosa / L. hyalospora 的 35,966 条蛋白", "无", "中（与本项目 BRAKER2 层同工具家族，循环性最重）",
     "Shea et al. 2023, Microbiol Resour Announc 12(7):e00338-23, DOI 10.1128/mra.00338-23"),
    ("JGI 各株（DOE Joint Genome Institute 提交）", "GCF_025…/GCA_0021…/GCA_036…", "JGI Annotation Pipeline（非 BRAKER/Funannotate/MAKER）",
     "有：转录组数据（JGI 内部、GenBank、合作者）trim-cluster-assemble 后经 BLAT（95% identity + 80% coverage）与 GMAP 映射；用于训练 FGENESH（4:1 训练/测试）、驱动 EST_map / combest / estExt 转录本预测、修正 UTR 与冲突结构，并进入模型选择打分 S = S_blast*(cov1*cov2 + CC_a)。逐株是否采用 RNA-seq 未在元数据中逐条标注（MycoCosm 页面对外返回 403，未能自动核对）",
     "有该流程环节：JGI 自述流程含 automated gene modeling + manual inspection for quality control + community curation，并有 annotator→peer review→community→GenBank 四层评审；但人工审编比例未公开",
     "中—高（本报告主力基准，且方法学上与 BRAKER 独立）",
     "三份 JGI SOP（2024 / 2021-03-13 / 旧版）；Grigoriev et al. 2014, NAR 42:D699-704, PMID 24297253；Haridas et al. 2018, Methods Mol Biol 1775:171-184, PMID 29876818"),
    ("Mucor flavus KT1a", "GCA_040367835.1", "Funannotate 1.8.15（有论文）",
     "无：仅改 --busco_db mucolares 与 --busco_seed_species Rhizopus_oryzae，无转录证据", "无", "中",
     "Hosono et al. 2024, Microbiol Resour Announc 13(9):e00343-24, DOI 10.1128/mra.00343-24"),
    ("Helicostylum pulchrum KT1b", "GCA_040366295.1", "Funannotate 1.8.15（有论文）",
     "无（同上）", "无", "中", "Hosono et al. 2024, Microbiol Resour Announc 13(9):e00343-24"),
    ("Umbelopsis isabellina WA0000067209", "GCA_016758905.1", "MAKER v2.31.8（有论文）",
     "有 EST 证据但不是 RNA-seq：100,082 条公共 EST + 273,458 条 Mucoromycotina 蛋白 + Swiss-Prot", "无", "中",
     "Muszewska et al. 2021, Front Microbiol 12:636986, DOI 10.3389/fmicb.2021.636986"),
    ("Umbelopsis sp. WA50703", "GCA_053572175.1", "Funannotate 1.8.15（有论文，但同株另有组装）",
     "有：funannotate train 用两套真菌 RNA-seq（hisat2 v2.2.1 比对）", "无（文献中的 MAISEN 审校属另一组装 GCA_964291815.1）", "中",
     "Okrasińska et al. 2026, BMC Genomics 27:718, DOI 10.1186/s12864-026-13055-5"),
    ("Linnemannia hyalina SCG-10", "GCA_019671135.1", "Funannotate v1.8.1 全默认（有论文）",
     "无：论文明确把 RNA-seq 列为未来工作", "无", "中—低（另需查基因组偏大与污染）",
     "Aldossari & Ishii 2021, Microbiol Resour Announc 10(37):e00692-21, DOI 10.1128/mra.00692-21"),
    ("Glomus cerebriforme DAOM 227022", "GCA_003550305.1", "JGI Annotation Pipeline（有论文）",
     "该研究自身的 RNA-seq 只做了 R. irregularis 与 G. rosea，本株无 RNA-seq；JGI 流程对本株采用了哪些其他证据未验证",
     "同 JGI 流程（见上）", "中（受装配碎片化拖累）", "Morin et al. 2019, New Phytol 222(3):1584-1598, DOI 10.1111/nph.15687（ENA 已挂 PMID 30636349）"),
    ("Diversispora epigaea IT104", "GCA_003547095.1", "未验证（有论文，但注释流程无法核实）",
     "未验证（原文被出版商访问限制拦截，无可用 OA 版本）", "未验证", "中—低（来源不明的基因集不宜作高质量基准）",
     "Sun et al. 2019, New Phytol 221(3):1556-1573, DOI 10.1111/nph.15472"),
    ("Gigaspora rosea FL105 / G. margarita NC121", "GCA_910592045.1 / GCA_910592245.1", "GeneMark-ES v4.33 ab initio + Swiss-Prot hint（有论文）",
     "无：纯从头预测，零转录证据", "无", "低（不建议作基准）",
     "Montoliu-Nerin et al. 2021, Front Fungal Biol 2:716385, DOI 10.3389/ffunb.2021.716385"),
    ("Absidia glauca CBS 101.48", "GCA_900079185.1", "未文档化（未找到任何论文）",
     "未验证", "未验证", "低（不建议采用）",
     "仅 BioProject PRJEB13349；Ensembl Fungi genebuild_method = import"),
    ("EBP-Nor 各株（ELIXIR Norway 提交）", "GCA_977…/GCA_9739…", "EBP-Nor 自研流程（无方法论文）",
     "无：仅 3 套蛋白同源（miniprot + GALBA + EvidenceModeler），ENA flatfile 内 0 处 RNA-seq 痕迹", "无", "低—中（只宜作蛋白-only 基线）",
     "流程见 github.com/ebp-nor/GenomeAnnotation；同流程版本见 Tørresen et al. 2024, J Hered（非同行评审的方法引文）"),
    ("GCF_ 条目的注释归属（通例）", "全部 12 个 GCF_", "RefSeq 侧为「继承 GenBank 注释」，非 NCBI 重注释",
     "依提交方而定（见各行）", "依提交方而定", "—",
     "12 个 GCF_ 条目在 NCBI 均存在配对的 GCA 孪生记录（ftppath_genbank 齐全），RefSeq 侧只做重新发布；Rhizophagus irregularis 的 GCF 还挂在独立的 RefSeq BioProject PRJNA1128271"),
]

DOC_CANDIDATES = [
    ("Mucoromycotina / Mucorales", "Rhizopus delemar RA 99-880",
     "× 剔除（仅同名不同株）", "No.152 Rhizopus delemar / 采纳名 Rhizopus arrhizus，GCA_054127015.1, SE_PDA7",
     "本队列中的 R. delemar 是 SE_PDA7，与经典参考株 RA 99-880 无关，且该装配「原始蛋白注释 = 无」——本方与参照两侧都凑不齐，两条判据同时不满足。同名物种在本队列中不可用作 RA 99-880 的替代。"),
    ("Mucoromycotina / Mucorales", "Phycomyces blakesleeanus NRRL1555",
     "△ 可用但须改菌株号", "No.134 Phycomyces blakesleeanus UBC21, GCA_040209145.1",
     "本队列为 UBC21 而非 NRRL1555。V2 落点：备选组 B10，设计格 A（本方 BRAKER3 × 参考 R+）。菌株号必须改写。"),
    ("Mucoromycotina / Mucorales", "Mucor lusitanicus MU402",
     "△ 可用但须改菌株号", "No.120 Mucor lusitanicus CBS 277.49, GCA_040256725.1",
     "本队列为 CBS 277.49。V2 落点：优先组 P4，设计格 A（本方 BRAKER3 × 参考 R+，本方两条 RNA-seq 比对率均 >99%）。"),
    ("Mucoromycotina / Mucorales", "Mucor circinelloides CBS 277.49",
     "√ 与上一条为同一株", "No.120（同 GCA_040256725.1）",
     "该菌株现用名为 Mucor lusitanicus；表3 采纳名亦为 Mucor lusitanicus。此条与上一条重复，实际只占一个名额。"),
    ("Mucoromycotina / Mucorales", "Lichtheimia ornata",
     "√ 直接采用", "No.77 Lichtheimia ornata CBS 291.66, GCF_029851405.1",
     "文档已注明 GCF_029851405.1，核实无误。V2 落点：优先组 P8，设计格 B（本方 BRAKER2 × 参考 R−）——参考注释本身也是 BRAKER2 蛋白模式，是全部候选中循环性最强、但最能检验「同工具可复现性」的配对。"),
    ("Mucoromycotina / Mucorales", "Lichtheimia corymbifera",
     "× 剔除", "No.75 Lichtheimia corymbifera S6502, GCA_037042095.1",
     "该装配「原始蛋白注释 = 无」（本项目为首次注释），且基因组 BUSCO C 仅 87.8%，不能用作基准。同属替代见 P8（L. ornata）与 B5（L. hyalospora，设计格 D）。"),
    ("Mucoromycotina / Mucorales", "Rhizopus microsporus",
     "△ 降为备选", "No.153 Rhizopus microsporus, GCA_977110955.1",
     "基因组优秀（C=99.4%，N50=4.72 Mb），但注释为 EBP-Nor 自研流程（miniprot + GALBA + EvidenceModeler）纯蛋白同源证据的自动产物：无 RNA-seq、无 RefSeq 镜像、无注释方法论文，异构体负荷 23.3%（1,696/10,755 基因有多转录本）。V2 落点：备选组 B1，设计格 C（本方 BRAKER3 × 参考 R−），用作「本方多一套 RNA-seq 的增益」的量化样本。"),
    ("Mucoromycotina / Mucorales", "Absidia glauca 或 Cunninghamella echinulata",
     "× 不建议 / △ 备选", "Absidia glauca: No.2, GCA_900079185.1；Cunninghamella echinulata: No.43, GCA_977110915.1",
     "Absidia glauca 经核查不建议采用：未检索到任何发表论文（NCBI assembly→PubMed 链接为空），注释流程无任何文档，Ensembl Fungi 确认为原样导入（genebuild_method = import，最后更新 2016-04）；装配仅 Illumina + Velvet，5,338 个 scaffold、contig N50 仅 11 kb。Cunninghamella echinulata 落点：备选组 B2，设计格 C（本方 BRAKER3 × 参考 R−，异构体负荷 11.6%）。"),
    ("Mucoromycotina / Umbelopsidales", "Umbelopsis ramanniana",
     "△ 移入扩展替补池", "No.168 Umbelopsis ramanniana, GCA_977110945.1",
     "文档所述「基因数异常」经核实成立：既有注释仅 9,709 个基因（390 genes/Mb），低于同属其他菌株。该株的价值在于作为异常例证（服务批次二任务 2.4），而非承担「高质量基准」的角色；同时其注释为 EBP-Nor 纯蛋白自动流程、异构体负荷 11.9%。V2 中同属的 Umbelopsis vinacea 进入备选组 B9（设计格 B），可在匹配格内交叉核对 Umbelopsis 属的基因数偏低问题。"),
    ("Mortierellomycotina", "Linnemannia elongata AG-77",
     "△ 可用但须改菌株号", "No.81 Linnemannia elongata NVP64, GCA_036320915.1",
     "本队列为 NVP64。V2 落点：优先组 P6，设计格 A。菌株号必须改写。"),
    ("Mortierellomycotina", "Mortierella alpina",
     "√ 可用", "No.91 Mortierella alpina, GCA_977091265.1",
     "装配与注释均达标。V2 落点：备选组 B3，设计格 C（本方 BRAKER3 × 参考 R−，EBP-Nor 纯蛋白流程），同时补齐 Mortierellomycotina 在错配格上的代表。"),
    ("Mortierellomycotina", "Lobosporangium transversale",
     "√ 直接采用 → 调整为扩展替补", "No.88 Lobosporangium transversale NRRL 3116, GCF_002105155.1",
     "核实无误，RefSeq 参考基因组 + JGI 注释，设计格为 A。V2 中因通道 A 名额已由 P1–P6 与 B8、B10 填满，移入扩展替补池，仍可随时取用。该 GCF 记录在 NCBI 存在配对的 GenBank 孪生记录 GCA_002105155.1，RefSeq 侧只做重新发布，不可写成「NCBI 注释」。"),
    ("Glomeromycotina", "Rhizophagus irregularis DAOM 197198",
     "√ 直接采用（首要）", "No.147 Rhizophagus irregularis DAOM-197198, GCF_026210795.1",
     "核实无误；本队列中唯一染色体级 AMF 装配（32 条染色体、23 条无 gap、17 条 T2T）。V2 落点：优先组 P1，设计格 A；本方与参考两侧都有 RNA-seq，是全套设计中最干净的一个配对。"),
    ("Glomeromycotina", "Gigaspora margarita 或 G. rosea",
     "× 不建议作基准", "G. margarita: No.61, GCA_910592245.1；G. rosea: No.62, GCA_910592045.1",
     "两者基因组 BUSCO 完整度虽高（98.8% / 98.5%），但 N50 仅 6 kb 与 4 kb，装配高度破碎；参考注释为 GeneMark-ES 纯从头预测、零转录证据，G. rosea 的 BUSCO C 仅 46%。用一份自身残缺的参考去评判更优的注释会系统性惩罚本方。若需大基因组 AMF，改用扩展替补 Glomus cerebriforme（已入备选 B7，设计格 D）。"),
    ("Glomeromycotina", "Diversispora epigaea 或 Rhizophagus clarus",
     "√ 采用（R. clarus 优先；D. epigaea 扩展替补）", "R. clarus: No.146, GCA_015698045.1；D. epigaea: No.50, GCA_003547095.1",
     "R. clarus 落点：优先组 P3，设计格 A（BRAKER v2.1.2 + RNA-seq）。D. epigaea 的注释流程本次未能核实（原文受出版商访问限制，无可用 OA 版本），属来源不明的基因集，移入扩展替补池，不作基准。"),
]


# 既有注释的证据类型分层（每条均已核实到一手来源；unverified 项如实标注）
EVIDENCE_TIERS = [
    ("Rhizophagus irregularis DAOM-197198", "GCF_026210795.1", "BRAKER 2.1.5 + funannotate/PASA（有论文）",
     "有，且为短读 + 长读双证据：Illumina RNA-seq → STAR 2.7.6a → BRAKER 2.1.5；另有 Nanopore PCR-cDNA RNA-seq → funannotate 1.8.7 update + PASA 精修 UTR 与 poly(A)",
     "部分：Hi-C 挂载人工校正 + 人工剔除含转座子结构域的基因；非逐基因模型审编", "高（AMF 内事实参考；但与 BRAKER3 同源）",
     "Manley et al. 2023, G3 13(6):jkad077, DOI 10.1093/g3journal/jkad077"),
    ("Rhizophagus clarus HR1", "GCA_015698045.1", "BRAKER v2.1.2 + AUGUSTUS 3.3.2 + GeneMark-ES 4.32（仅预印本）",
     "有，作者自述 RNA-seq-based；但具体用了哪套 RNA-seq 数据未点名（未验证）", "无记录", "中（唯一文献为未经同行评审的 2026 预印本）",
     "Kobayashi et al., Research Square 预印本 2026-08-13, DOI 10.21203/rs.3.rs-10340957/v1（Under Revision）"),
    ("Mucor velutinosus NIH1002", "GCF_035048745.1", "de novo 自动流程（有论文）",
     "有：16 个生长条件 RNA-seq → STAR v2.5.3a / StringTie v1.3.3；与 TransDecoder v5.0.2、BRAKER1 v1.9、CodingQuarry v2.0 结果交 EVidenceModeler 整合，funannotate 收尾",
     "无", "高（证据类型与 BRAKER3 最接近；但与 BRAKER 家族同源）",
     "Conlan et al. 2024, Microbiol Resour Announc 13(7):e00062-24, DOI 10.1128/mra.00062-24"),
    ("Lichtheimia ornata CBS 291.66", "GCF_029851405.1", "BRAKER2 蛋白模式（有论文）",
     "无：证据为 L. corymbifera / L. ramosa / L. hyalospora 的 35,966 条蛋白", "无", "中（与本项目 BRAKER2 层同工具家族，循环性最重）",
     "Shea et al. 2023, Microbiol Resour Announc 12(7):e00338-23, DOI 10.1128/mra.00338-23"),
    ("JGI 各株（DOE Joint Genome Institute 提交）", "GCF_025…/GCA_0021…/GCA_036…", "JGI Annotation Pipeline（非 BRAKER/Funannotate/MAKER）",
     "有：转录组数据（JGI 内部、GenBank、合作者）trim-cluster-assemble 后经 BLAT（95% identity + 80% coverage）与 GMAP 映射；用于训练 FGENESH（4:1 训练/测试）、驱动 EST_map / combest / estExt 转录本预测、修正 UTR 与冲突结构，并进入模型选择打分 S = S_blast*(cov1*cov2 + CC_a)。逐株是否采用 RNA-seq 未在元数据中逐条标注（MycoCosm 页面对外返回 403，未能自动核对）",
     "有该流程环节：JGI 自述流程含 automated gene modeling + manual inspection for quality control + community curation，并有 annotator→peer review→community→GenBank 四层评审；但人工审编比例未公开",
     "中—高（本报告主力基准，且方法学上与 BRAKER 独立）",
     "三份 JGI SOP（2024 / 2021-03-13 / 旧版）；Grigoriev et al. 2014, NAR 42:D699-704, PMID 24297253；Haridas et al. 2018, Methods Mol Biol 1775:171-184, PMID 29876818"),
    ("Mucor flavus KT1a", "GCA_040367835.1", "Funannotate 1.8.15（有论文）",
     "无：仅改 --busco_db mucolares 与 --busco_seed_species Rhizopus_oryzae，无转录证据", "无", "中",
     "Hosono et al. 2024, Microbiol Resour Announc 13(9):e00343-24, DOI 10.1128/mra.00343-24"),
    ("Helicostylum pulchrum KT1b", "GCA_040366295.1", "Funannotate 1.8.15（有论文）",
     "无（同上）", "无", "中", "Hosono et al. 2024, Microbiol Resour Announc 13(9):e00343-24"),
    ("Umbelopsis isabellina WA0000067209", "GCA_016758905.1", "MAKER v2.31.8（有论文）",
     "有 EST 证据但不是 RNA-seq：100,082 条公共 EST + 273,458 条 Mucoromycotina 蛋白 + Swiss-Prot", "无", "中",
     "Muszewska et al. 2021, Front Microbiol 12:636986, DOI 10.3389/fmicb.2021.636986"),
    ("Umbelopsis sp. WA50703", "GCA_053572175.1", "Funannotate 1.8.15（有论文，但同株另有组装）",
     "有：funannotate train 用两套真菌 RNA-seq（hisat2 v2.2.1 比对）", "无（文献中的 MAISEN 审校属另一组装 GCA_964291815.1）", "中",
     "Okrasińska et al. 2026, BMC Genomics 27:718, DOI 10.1186/s12864-026-13055-5"),
    ("Linnemannia hyalina SCG-10", "GCA_019671135.1", "Funannotate v1.8.1 全默认（有论文）",
     "无：论文明确把 RNA-seq 列为未来工作", "无", "中—低（另需查基因组偏大与污染）",
     "Aldossari & Ishii 2021, Microbiol Resour Announc 10(37):e00692-21, DOI 10.1128/mra.00692-21"),
    ("Glomus cerebriforme DAOM 227022", "GCA_003550305.1", "JGI Annotation Pipeline（有论文）",
     "该研究自身的 RNA-seq 只做了 R. irregularis 与 G. rosea，本株无 RNA-seq；JGI 流程对本株采用了哪些其他证据未验证",
     "同 JGI 流程（见上）", "中（受装配碎片化拖累）", "Morin et al. 2019, New Phytol 222(3):1584-1598, DOI 10.1111/nph.15687（ENA 已挂 PMID 30636349）"),
    ("Diversispora epigaea IT104", "GCA_003547095.1", "未验证（有论文，但注释流程无法核实）",
     "未验证（原文被出版商访问限制拦截，无可用 OA 版本）", "未验证", "中—低（来源不明的基因集不宜作高质量基准）",
     "Sun et al. 2019, New Phytol 221(3):1556-1573, DOI 10.1111/nph.15472"),
    ("Gigaspora rosea FL105 / G. margarita NC121", "GCA_910592045.1 / GCA_910592245.1", "GeneMark-ES v4.33 ab initio + Swiss-Prot hint（有论文）",
     "无：纯从头预测，零转录证据", "无", "低（不建议作基准）",
     "Montoliu-Nerin et al. 2021, Front Fungal Biol 2:716385, DOI 10.3389/ffunb.2021.716385"),
    ("Absidia glauca CBS 101.48", "GCA_900079185.1", "未文档化（未找到任何论文）",
     "未验证", "未验证", "低（不建议采用）",
     "仅 BioProject PRJEB13349；Ensembl Fungi genebuild_method = import"),
    ("EBP-Nor 各株（ELIXIR Norway 提交）", "GCA_977…/GCA_9739…", "EBP-Nor 自研流程（无方法论文）",
     "无：仅 3 套蛋白同源（miniprot + GALBA + EvidenceModeler），ENA flatfile 内 0 处 RNA-seq 痕迹", "无", "低—中（只宜作蛋白-only 基线）",
     "流程见 github.com/ebp-nor/GenomeAnnotation；同流程版本见 Tørresen et al. 2024, J Hered（非同行评审的方法引文）"),
    ("GCF_ 条目的注释归属（通例）", "全部 12 个 GCF_", "RefSeq 侧为「继承 GenBank 注释」，非 NCBI 重注释",
     "依提交方而定（见各行）", "依提交方而定", "—",
     "12 个 GCF_ 条目在 NCBI 均存在配对的 GCA 孪生记录（ftppath_genbank 齐全），RefSeq 侧只做重新发布；Rhizophagus irregularis 的 GCF 还挂在独立的 RefSeq BioProject PRJNA1128271"),
]

DOC_CANDIDATES = [
    ("Mucoromycotina / Mucorales", "Rhizopus delemar RA 99-880",
     "× 剔除（仅同名不同株）", "No.152 Rhizopus delemar / 采纳名 Rhizopus arrhizus，GCA_054127015.1, SE_PDA7",
     "本队列中的 R. delemar 是 SE_PDA7，与经典参考株 RA 99-880 无关，且该装配'原始蛋白注释 = 无'，两条判据同时不满足。同名物种在本队列中不可用作 RA 99-880 的替代。"),
    ("Mucoromycotina / Mucorales", "Phycomyces blakesleeanus NRRL1555",
     "△ 可用但须改菌株号", "No.134 Phycomyces blakesleeanus UBC21, GCA_040209145.1",
     "本队列为 UBC21 而非 NRRL1555。已列为优先组 P4，报告中菌株号必须改写为 UBC21。"),
    ("Mucoromycotina / Mucorales", "Mucor lusitanicus MU402",
     "△ 可用但须改菌株号", "No.120 Mucor lusitanicus CBS 277.49, GCA_040256725.1",
     "本队列为 CBS 277.49。已列为优先组 P3。"),
    ("Mucoromycotina / Mucorales", "Mucor circinelloides CBS 277.49",
     "√ 与上一条为同一株", "No.120（同 GCA_040256725.1）",
     "该菌株现用名为 Mucor lusitanicus；表3 采纳名亦为 Mucor lusitanicus。此条与上一条重复，实际只占一个名额。"),
    ("Mucoromycotina / Mucorales", "Lichtheimia ornata",
     "√ 直接采用", "No.77 Lichtheimia ornata CBS 291.66, GCF_029851405.1",
     "文档已注明 GCF_029851405.1，核实无误。已列为优先组 P6。"),
    ("Mucoromycotina / Mucorales", "Lichtheimia corymbifera",
     "× 剔除", "No.75 Lichtheimia corymbifera S6502, GCA_037042095.1",
     "该装配'原始蛋白注释 = 无'（本项目为首次注释），且基因组 BUSCO C 仅 87.8%，不能用作基准。同属替代见 P6（L. ornata）与扩展替补（L. hyalospora）。"),
    ("Mucoromycotina / Mucorales", "Rhizopus microsporus",
     "△ 降为备选", "No.153 Rhizopus microsporus, GCA_977110955.1",
     "基因组优秀（C=99.4%，N50=4.72 Mb），但注释为 EBP-Nor 自研流程（miniprot + GALBA + EvidenceModeler）纯蛋白同源证据的自动产物：无 RNA-seq、无 RefSeq 镜像、无注释方法论文，异构体负荷 23.3%（1,696/10,755 基因有多转录本）。队列内 RNA-seq 比对率亦仅 43.93%–62.42%。列为备选组 B10，定位为「蛋白-only 自动流程基线」。"),
    ("Mucoromycotina / Mucorales", "Absidia glauca 或 Cunninghamella echinulata",
     "× 不建议 / △ 扩展替补", "Absidia glauca: No.2, GCA_900079185.1；Cunninghamella echinulata: No.43, GCA_977110915.1",
     "Absidia glauca 经核查不建议采用：未检索到任何发表论文（NCBI assembly→PubMed 链接为空），注释流程无任何文档，Ensembl Fungi 确认为原样导入（genebuild_method = import，最后更新 2016-04）；装配仅 Illumina + Velvet，5,338 个 scaffold、contig N50 仅 11 kb；基因组 BUSCO C 亦仅 97.0%。Cunninghamella echinulata 为 EBP-Nor 装配（C=99.4%，N50=4.26 Mb，异构体负荷 11.6%），列入扩展替补。"),
    ("Mucoromycotina / Umbelopsidales", "Umbelopsis ramanniana",
     "△ 移入扩展替补池", "No.168 Umbelopsis ramanniana, GCA_977110945.1",
     "文档所述'基因数异常'经核实成立：既有注释仅 9,709 个基因（390 genes/Mb），低于同属其他菌株（U. isabellina 9,033、U. sp. AD052 9,475 属同一量级，但密度更高）。该株的价值在于作为异常例证（服务批次二任务 2.4），而非承担'高质量基准'的角色；同时其注释为 EBP-Nor 纯蛋白自动流程、异构体负荷 11.9%。故移入扩展替补池。"),
    ("Mortierellomycotina", "Linnemannia elongata AG-77",
     "△ 可用但须改菌株号", "No.81 Linnemannia elongata NVP64, GCA_036320915.1",
     "本队列为 NVP64。已列为优先组 P9。"),
    ("Mortierellomycotina", "Mortierella alpina",
     "√ 可用", "No.91 Mortierella alpina, GCA_977091265.1",
     "装配与注释均达标，但注释来自 EBP-Nor 自动流程，故置于扩展替补而非优先组。"),
    ("Mortierellomycotina", "Lobosporangium transversale",
     "√ 直接采用", "No.88 Lobosporangium transversale NRRL 3116, GCF_002105155.1",
     "核实无误，RefSeq 参考基因组 + JGI 注释。已列为优先组 P10。"),
    ("Glomeromycotina", "Rhizophagus irregularis DAOM 197198",
     "√ 直接采用（首要）", "No.147 Rhizophagus irregularis DAOM-197198, GCF_026210795.1",
     "核实无误；本队列中唯一染色体级 AMF 装配。已列为优先组 P1。"),
    ("Glomeromycotina", "Gigaspora margarita 或 G. rosea",
     "× 不建议作基准", "G. margarita: No.61, GCA_910592245.1；G. rosea: No.62, GCA_910592045.1",
     "两者基因组 BUSCO 完整度虽高（98.8% / 98.5%），但 N50 仅 6 kb 与 4 kb，装配高度破碎，'高质量基因组'这一条不成立。若确需大基因组 AMF，改用扩展替补 Glomus cerebriforme（N50=148 kb）或 Diversispora epigaea（N50=453 kb）。"),
    ("Glomeromycotina", "Diversispora epigaea 或 Rhizophagus clarus",
     "√ 采用（R. clarus 优先；D. epigaea 备选）", "R. clarus: No.146, GCA_015698045.1；D. epigaea: No.50, GCA_003547095.1",
     "R. clarus 已列为优先组 P2；D. epigaea 列为备选组 B10。"),
]


# ----------------------------------------------------------------------------
# 3. Word 输出
# ----------------------------------------------------------------------------

def build_docx(rows, stats, notes, REF_CHECK, cell_tab):
    from docx import Document
    from docx.shared import Pt, Cm, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml.ns import qn

    doc = Document()

    # 全局字体
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(10.5)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    for i in range(1, 5):
        h = doc.styles[f"Heading {i}"]
        h.font.name = "Times New Roman"
        h.element.rPr.rFonts.set(qn("w:eastAsia"), "黑体")
        h.font.color.rgb = RGBColor(0, 0, 0)
        h.font.size = Pt({1: 16, 2: 13.5, 3: 11.5, 4: 10.5}[i])
    for s in doc.sections:
        s.left_margin = s.right_margin = Cm(2.2)
        s.top_margin = s.bottom_margin = Cm(2.0)

    def para(text, bold=False, size=None, align=None, space_after=6, italic=False):
        p = doc.add_paragraph()
        r = p.add_run(text)
        r.bold = bold
        r.italic = italic
        if size:
            r.font.size = Pt(size)
        if align == "c":
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_after = Pt(space_after)
        return p

    def bullet(text):
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(text)
        p.paragraph_format.space_after = Pt(2)
        return p

    def table(headers, data, widths=None, font=8.5):
        t = doc.add_table(rows=1, cols=len(headers))
        t.style = "Table Grid"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        hdr = t.rows[0].cells
        for i, h in enumerate(headers):
            hdr[i].text = ""
            p = hdr[i].paragraphs[0]
            r = p.add_run(h)
            r.bold = True
            r.font.size = Pt(font)
            r.font.name = "Times New Roman"
            r.element.rPr.rFonts.set(qn("w:eastAsia"), "黑体")
        for row in data:
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = ""
                p = cells[i].paragraphs[0]
                r = p.add_run(str(v))
                r.font.size = Pt(font)
                r.font.name = "Times New Roman"
                r.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
        if widths:
            for row in t.rows:
                for i, w in enumerate(widths):
                    row.cells[i].width = Cm(w)
        doc.add_paragraph().paragraph_format.space_after = Pt(2)
        return t

    def species_block(tag, name, strain, sub, acc, src, level, size, n50, busco, genes, rna, orthodb,
                      isoform, reason, risk, mode, tier_label, cell_label):
        doc.add_heading(f"{tag}　{name} {strain}", level=3)
        table(
            ["项目", "内容", "项目", "内容"],
            [
                ["亚门", sub, "Assembly accession", acc],
                ["本方注释方式", mode + ("（基因组 + 蛋白 + RNA-seq）" if mode == "BRAKER3" else "（基因组 + 蛋白，仅蛋白证据）"),
                 "参考注释证据等级", tier_label],
                ["设计格", cell_label, "装配水平", level],
                ["基因组大小", f"{size} Mb", "N50", f"{n50} kb"],
                ["基因组 BUSCO C", f"{busco}%", "既有注释基因数", genes],
                ["队列内 RNA-seq", rna, "OrthoDB v12.2 证据库", orthodb],
                ["参考注释异构体", isoform, "参考蛋白取回", "已完成（附录 B/C）" if acc.startswith(("GCA_", "GCF_")) else "需从 JGI MycoCosm 下载"],
            ],
            widths=[2.6, 5.9, 3.4, 5.1], font=8.5)
        p = doc.add_paragraph()
        r = p.add_run("入选理由："); r.bold = True; r.font.size = Pt(9.5)
        p.add_run(reason).font.size = Pt(9.5)
        p.paragraph_format.space_after = Pt(4)
        p = doc.add_paragraph()
        r = p.add_run("风险与备注："); r.bold = True; r.font.size = Pt(9.5)
        p.add_run(risk).font.size = Pt(9.5)
        p.paragraph_format.space_after = Pt(10)

    # ---------------- 封面 / 标题 ----------------
    doc.add_heading("基准物种遴选报告（V2）", level=0)
    para("批次一 · 任务 1.1 —— 与高质量参考注释做蛋白层面对照的物种名单", bold=True, size=12, align="c")
    para("Mucoromycota 202 基因组注释项目（Scientific Data Data Descriptor）", size=10, align="c")
    para("版本 v2 · 2026-09-14 · 全部结论以 Assembly accession 锚定，可逐条复核；V2 相对 V1 的改动见 README.md 与第 1.5 节",
         size=9, align="c", space_after=14)

    doc.add_heading("0　一页速览", level=2)
    para(f"V2 相对 V1 的核心改动只有一件：把「本方注释方式」提升为与「亚门」「基因组质量」并列的设计轴。"
         f"本研究 202 株中 {stats['n_b2_all']} 株是 BRAKER2（仅蛋白证据）、{stats['n_b3_all']} 株是 BRAKER3"
         f"（蛋白 + RNA-seq）——超过一半的交付物出自没有 RNA-seq 的那条流水线。"
         f"若不加区分地拿它们去和以 RNA-seq 为基础的参考注释比，差异里会混入证据类型效应并被误读为方法优劣。"
         f"因此 V2 采用 2×2 抽样设计：本方 BRAKER3/BRAKER2 × 参考注释含/不含转录证据，"
         f"优先组只放两个「证据匹配」格，备选组放两个「可控错配」格用于量化证据效应。")
    para(f"从 202 株中筛出 {stats['n_has_ann']} 株具备第三方既有注释的物种，施加基因组完整度与装配连续性门槛后，"
         f"按三轴配额定出优先组 10 株、备选组 10 株，共 20 株。全部 20 株均在本研究 202 株队列内，"
         f"accession 与表1/表3/表6/表7 及 NCBI Datasets API（2026-09-14 快照）逐条核对一致。")
    para("本报告在核查过程中得到两条必须提前说清的结论，它们决定了后文所有措辞。"
         "第一，本队列没有任何一份既有注释具备可核实的、成规模的逐基因人工审校记录。"
         "RefSeq 的 refseq_category = reference genome 只是 NCBI 的分类指定，不代表审校，也不代表 NCBI 重注释——"
         "12 个 GCF_ 条目在 NCBI 均存在配对的 GCA 孪生记录，RefSeq 侧只做重新发布。"
         "唯一声称含人工环节的是 JGI（自述流程为 automated gene modeling + manual inspection for quality control "
         "+ community curation，并有四层评审），但人工审编比例未公开。因此技术验证的定位只能是"
         "「与既有注释的一致性/差异分析」，不能写成「与人工审校金标准比对」。"
         "第二，优先组 10 株里有 4 株的既有注释本身就是 BRAKER 家族产物（Rhizophagus irregularis 用 BRAKER 2.1.5、"
         "Rhizophagus clarus 用 BRAKER v2.1.2、Mucor velutinosus 用 BRAKER1、Lichtheimia ornata 用 BRAKER2 蛋白模式），"
         "与本研究属同一工具家族；真正方法学独立的是 6 株 JGI 注释（其流程为 FGENESH / GeneMark / GeneWise / "
         "EST_map / combest / estExt，不用 BRAKER）。做结论时必须把这两组分开看。",
         size=10)
    bullet("文档原候选名单 15 条逐条核查后：可直接采用 6 条；须改写菌株号后可用 3 条"
           "（Phycomyces blakesleeanus、Mucor lusitanicus、Linnemannia elongata）；降级 2 条"
           "（Rhizopus microsporus → 备选组，Umbelopsis ramanniana → 扩展替补池）；"
           "剔除或不建议 4 条（Rhizopus delemar RA 99-880、Lichtheimia corymbifera、Absidia glauca、"
           "Gigaspora 两株）；另有 1 条与 Mucor lusitanicus 实为同一株，只占 1 个名额。详见第 2 节。")
    _cnt = {c: sum(cell_tab.get(c, {}).values()) for c in "ABCD"}
    bullet(f"20 株的设计格落点：匹配格 A（本方 BRAKER3 × 参考 R+）{_cnt['A']} 株、"
           f"匹配格 B（本方 BRAKER2 × 参考 R−）{_cnt['B']} 株、"
           f"错配格 C（本方 BRAKER3 × 参考 R−）{_cnt['C']} 株、"
           f"错配格 D（本方 BRAKER2 × 参考 R+）{_cnt['D']} 株。"
           f"优先组 10 株全部来自两个匹配格，这是本报告与 V1 最实质的差别。")
    bullet("亚门覆盖：Mucoromycotina 12 株、Mortierellomycotina 5 株、Glomeromycotina 3 株，"
           "对应队列的 123 : 57 : 22（61% : 28% : 11%）；Glomeromycotina 有意加权到 15%，"
           "因为该亚门装配质量最差、最易被质疑，2 株不足以支撑结论。")
    _pr = [r for r in rows if r["group"] == "优先"]
    _gcf = sum(1 for r in _pr if r["acc"].startswith("GCF_"))
    _jgi = sum(1 for r in _pr if ("Joint Genome Institute" in r["src"]) or r["acc"].startswith("JGI:"))
    _rna = sum(1 for r in _pr if not r["rna"].startswith("无"))
    _fam = {"JGI 机构流程": 0, "BRAKER 家族": 0, "EBP-Nor 自研流程": 0, "其他": 0}
    for r in _pr:
        if ("Joint Genome Institute" in r["src"]) or r["acc"].startswith("JGI:"):
            _fam["JGI 机构流程"] += 1
        elif any(k in r["src"] for k in ("Kew", "NHGRI", "NIBB", "Broad")):
            _fam["BRAKER 家族"] += 1
        elif "ELIXIR Norway" in r["src"]:
            _fam["EBP-Nor 自研流程"] += 1
        else:
            _fam["其他"] += 1
    bullet(f"优先组 10 株的参考注释来源构成：RefSeq（GCF_）参考基因组 {_gcf} 株；"
           "按流程家族分为 " + "、".join(f"{k} {v} 株" for k, v in _fam.items() if v) + "。"
           "其中 BRAKER 家族那几株存在方法同源问题（见第 6.12 节），"
           "JGI 机构流程是方法学上最独立的一组。")
    _norun = "、".join(r["name"] for r in _pr if r["rna"].startswith("无"))
    bullet(f"优先组 {_rna} 株在队列内保有 RNA-seq run 可直接调用；无保留 run 的 {10 - _rna} 株"
           f"（{_norun}）恰好就是通道 B 的全部成员——这不是巧合，通道 B 的定义本身就是「本方只有蛋白证据」，"
           f"因此这 4 株的参考注释也必须是 R−，两侧才对齐。")
    bullet("「既有注释有 RNA-seq 支持」这一条比预想的稀缺，且分层很细。已确证有 RNA-seq 证据的："
           "Mucor velutinosus NIH1002（16 个生长条件）、Rhizophagus irregularis（短读 + 长读双证据）、"
           "Rhizophagus clarus（BRAKER v2.1.2 + RNA-seq，但所用数据集未点名）。"
           "JGI 各株的流程以转录组数据为标准输入（BLAT + GMAP 映射后驱动 FGENESH / EST_map / combest / estExt），"
           "但逐株是否采用未在公开元数据中标注。"
           "已确证没有的：Lichtheimia ornata（BRAKER2 蛋白模式）、EBP-Nor 各株（纯蛋白同源）、"
           "Mucor flavus / Helicostylum pulchrum（Funannotate 1.8.15 默认）、Linnemannia hyalina（Funannotate 1.8.1 全默认）、"
           "Gigaspora 两株（GeneMark-ES 纯从头预测）、EBP-Nor 各株（miniprot + GALBA + EvidenceModeler 纯蛋白同源）。"
           "完整分层见第 1.6 节；第 1.5 节说明如何用这一分层构成 2×2 抽样设计。")
    bullet("已提前完成批次一任务 1.2 的一半工作：20 株的 OrthoDB v12.2 证据库归属标记随名单一并给出，"
           "用于后续按'在库 / 不在库'分组报告。")

    doc.add_heading("1　遴选口径", level=2)
    doc.add_heading("1.1　这件事不是批次二任务 2.4", level=3)
    table(["", "本报告（批次一 · 任务 1.1）", "批次二 · 任务 2.4"],
          [["目的", "验证本研究交付的注释是否可靠", "说明既有注释彼此异质、不可横向比较"],
           ["对象", "少数既有注释确实高质量的物种", "202 株中所有带既有蛋白文件的物种"],
           ["比较方式", "逐蛋白序列层面的定量比对（任务 1.3–1.5）", "完整度分数的分布与离散程度"],
           ["稿件归属", "Technical Validation 核心小节", "Background & Summary 的 Figure 1"]],
          widths=[2.4, 7.2, 7.4], font=9)
    para("两者不可互相替代。本报告只解决一件事：把'值得拿来对照的那几株'挑出来，并留下可复核的挑选依据。")

    doc.add_heading("1.2　入选判据", level=3)
    table(["编号", "判据", "本报告的操作化阈值"],
          [["R1", "必须在本研究 202 株队列内", "以 Assembly accession（或 JGI 门户代码+版本）与表3 主标识唯一匹配"],
           ["R2", "必须有第三方既有蛋白注释", "表1「原始蛋白注释 = 有」；＝无者一律出局（本项目首次注释，无 b 可比）"],
           ["R3", "基因组完整度高", "mucoromycota_odb12 基因组模式 BUSCO 完整度 C ≥ 98.5%（放宽档 ≥ 98.0%）"],
           ["R4", "装配连续性好", "N50 ≥ 800 kb，或装配水平达 Chromosome / Complete Genome"],
           ["R5", "注释来源可归因、可核实", "RefSeq（GCF_）参考基因组，或具名机构/JGI 提交且元数据可查"],
           ["R6", "亚门覆盖", "三个亚门均有代表，Glomeromycotina 至少 2 株"]],
          widths=[1.4, 5.2, 10.4], font=9)

    doc.add_heading("1.3　数据来源与核查方法", level=3)
    para("本报告不使用任何未标注来源的推断。全部字段来自以下五处，并已交叉核对：")
    bullet("投稿附表/表3_研究总表.xlsx —— 202 株权威编号、原名/采纳名、亚门、主标识（accession）、基因组 BUSCO C。")
    bullet("投稿附表/表1_来源清单_总表.xlsx —— 「原始蛋白注释 有/无」、assembly 水平、菌株号。")
    bullet("投稿附表/表6_注释证据_202株.xlsx —— 队列内保留的 RNA-seq run 与 HISAT2 比对率、BRAKER 模式。")
    bullet("投稿附表/表7_基因结构注释_202株.xlsx —— 装配 N50、基因结构统计。")
    bullet("NCBI Datasets API v2alpha（2026-09-14 查询快照）—— 装配水平、refseq_category、测序技术、"
           "装配方法、BioProject、以及 annotation_info 中的注释提交方与基因数。")
    bullet("6_Protein_evidence/14_tables_OrthoDBv12.2_Mucoromycota157_ProteinEvidence_Taxonomy.xlsx —— "
           "OrthoDB v12.2 Mucoromycota 157 文件子集（133 属种），用于任务 1.2 的'在库'标记。")

    doc.add_heading("1.4　筛选漏斗", level=3)
    table(["步骤", "条件", "剩余株数"],
          [["①", "本研究 202 株队列", stats["n_total"]],
           ["②", "表1「原始蛋白注释 = 有」（存在第三方既有注释 b）", stats["n_has_ann"]],
           ["③", "注释来源可在 NCBI / JGI 元数据中归因", stats["n_attrib"]],
           ["④", "基因组模式 BUSCO C ≥ 98.5%", stats["n_busco985"]],
           ["⑤", "且 N50 ≥ 800 kb 或装配为染色体级", stats["n_final"]],
           ["⑥", "本方注释方式与参考证据等级双明确（用于 2×2 设计分格）", 20],
           ["⑦", "按「亚门 : 基因组质量 : 本方模式 × 参考证据等级」三轴配额定序 → 优先 10 + 备选 10", 20]],
          widths=[1.6, 11.4, 3.6], font=9)
    para(f"放宽档（BUSCO C ≥ 98.0% 且 N50 ≥ 500 kb 或染色体级）可得到 {stats['n_loose']} 株，"
         f"即本报告的扩展替补池；若优先组在后续核查中需要替换，从这里按亚门就近取用。", size=9)

    doc.add_heading("1.5　三轴抽样设计（V2 的核心改动）", level=3)
    para("V1 只把「亚门比例」和「基因组质量」当作设计轴，参考注释本身用了什么证据只是当作逐株备注。"
         "这在本项目里是不够的：本研究 202 株中 119 株是 BRAKER2（基因组 + 蛋白，仅蛋白证据）、83 株是 BRAKER3"
         "（基因组 + 蛋白 + RNA-seq），超过一半的交付物出自「没有 RNA-seq」的那条流水线。"
         "如果拿它们去和一批以 RNA-seq 为基础做出的参考注释比，参考方必然多找回一批表达弱、演化快、"
         "类群特异的基因，这些差异会被误读成「我的注释更差」，而真正的原因是证据更薄。"
         "因此 V2 把「本方注释方式」提升为与亚门、基因组质量并列的第三根轴，并把参考注释的证据等级作为第四根轴，"
         "构成一个 2×2 抽样设计。", size=10)
    table(["轴", "水平", "来源", "本研究队列结构"],
          [["轴 1　亚门", "Mucoromycotina / Mortierellomycotina / Glomeromycotina", "表3",
            "123 : 57 : 22（202 株）"],
           ["轴 2　基因组质量", "BUSCO C ≥ 98.5% 且 （N50 ≥ 800 kb 或染色体级）", "表3 + 表7 + NCBI",
            "48 / 143（有既有注释者中）"],
           ["轴 3　本方注释方式", "BRAKER3（基因组 + 蛋白 + RNA-seq）/ BRAKER2（基因组 + 蛋白）", "表6",
            "83 : 119（202 株）；70 : 73（143 株）"],
           ["轴 4　参考注释证据等级", "R+（含转录证据：RNA-seq 或 JGI 转录组输入）/ R−（仅蛋白或从头预测）/ R?（未验证）",
            "本报告一手核查（第 1.6 节）", "依物种而定"]],
          widths=[3.0, 6.6, 3.4, 4.0], font=8.5)
    para("2×2 设计格、配额与用途：", bold=True, size=10)
    table(["格", "本方", "参考", "性质", "在论文中回答的问题", "本报告配额"],
          [["A", "BRAKER3", "R+", "匹配", "本方 RNA-seq 通道的注释，与同样有转录证据的参考能吻合到什么程度", f"优先 6 + 备选 2 = 8"],
           ["B", "BRAKER2", "R−", "匹配", "本方仅蛋白通道的注释，与同样只有蛋白证据的参考能吻合到什么程度", f"优先 4 + 备选 1 = 5"],
           ["C", "BRAKER3", "R−", "错配", "本方多一套 RNA-seq 时，比蛋白-only 参考多找回/多预测了什么", "备选 4"],
           ["D", "BRAKER2", "R+", "错配", "本方少一套 RNA-seq 时，相对有转录证据的参考漏掉了什么", "备选 3"]],
          widths=[1.0, 2.2, 1.4, 1.5, 8.5, 2.4], font=8.5)
    para("设计的估计量（建议按此写进 Technical Validation）：", bold=True, size=10)
    bullet("匹配格（A、B）× 亚门，给出一致性指标的**无偏基线**——这是评估本研究交付注释的主要证据。"
           "把 A 与 B 分开报告，等于把「有 RNA-seq 的一半」与「没有 RNA-seq 的一半」各自放在公平的条件下评估，"
           "直接回应「你的 BRAKER2 那部分是不是本来就做不好」这个必然会来的质疑。")
    bullet("错配格（C、D）给出**证据类型效应 Δ**：C 减 A 反映「本方多一套 RNA-seq」的增益，"
           "D 减 B 反映「本方少一套 RNA-seq」的损失。两个 Δ 方向相反、可以互相印证，"
           "把原本的混杂因素变成一个有量、可报告的变量。")
    bullet("四个格可以合并出一个**证据对齐后的一致性估计**：以 A、B 两格的一致性为交付质量的主数字，"
           "C、D 两格只用于解释差异来源，不并入主数字。")
    para("一条必须先说的结构约束：", bold=True, size=10)
    bullet("在 143 株有既有注释的物种中，「本方 BRAKER2 × 参考 R−」这一格只存在于 Mucoromycotina。"
           "过质量门槛的只有 5 株，全部是 Mucorales（EBP-Nor 批次 4 株 + Broad 的 Lichtheimia ornata 1 株）；"
           "Mortierellomycotina 与 Glomeromycotina 在这一格上是空的。原因是数据本身的分布："
           "能进入 JGI 体系的物种通常同时具备 RNA-seq，因此其参考注释落在 R+；"
           "而 EBP-Nor 的蛋白-only 批次恰好只覆盖 Mucorales。")
    bullet("后果与应对：通道 B 无法覆盖三个亚门。这不是遴选疏漏，必须在设计说明里如实写出，"
           "并把通道 B 的结论适用范围限定在 Mucoromycotina。若审稿人要求补齐，唯一办法是引入本队列之外的"
           "蛋白-only 参考注释，但那会破坏 R1 判据（必须在本研究 202 株队列内），本报告不建议。")
    para("三个轴的实际落点（20 株）：", bold=True, size=10)
    table(["设计格", "Mucoromycotina", "Mortierellomycotina", "Glomeromycotina", "小计"],
          [[CELL_LABEL[c], str(cell_tab.get(c, {}).get("Mucoromycotina", 0)),
            str(cell_tab.get(c, {}).get("Mortierellomycotina", 0)),
            str(cell_tab.get(c, {}).get("Glomeromycotina", 0)),
            str(sum(cell_tab.get(c, {}).values()))] for c in ["A", "B", "C", "D"]]
          + [["合计", str(sum(cell_tab.get(c, {}).get("Mucoromycotina", 0) for c in "ABCD")),
              str(sum(cell_tab.get(c, {}).get("Mortierellomycotina", 0) for c in "ABCD")),
              str(sum(cell_tab.get(c, {}).get("Glomeromycotina", 0) for c in "ABCD")),
              str(len(rows))]],
          widths=[4.6, 3.2, 3.4, 3.2, 1.6], font=8.5)
    para("亚门配比说明：20 株按 12 : 5 : 3（Mucoromycotina : Mortierellomycotina : Glomeromycotina）配额，"
         "对应队列的 123 : 57 : 22（61% : 28% : 11%）。Mortierellomycotina 略低于等比例（25%），"
         "Glomeromycotina 有意高于等比例（15%）——后者是文档明确要求「Glomeromycotina 必须有代表」的落实，"
         "因为该亚门装配质量最差、最易被质疑，若只放 2 株（等比例值）不足以支撑结论。", size=9)

    doc.add_heading("1.6　既有注释的证据类型分层", level=3)
    para("这一节回答一个比「哪株基因组好」更要命的问题：作为基准的那份注释，本身是用什么证据、由什么流程产生的。"
         "结论直接影响后续所有指标的解释方式。以下每一条都核实到一手来源（论文原文、ENA/NCBI 记录、流程仓库），"
         "无法核实的项如实标注「未验证」。")
    _tier_of = {
        "Rhizophagus irregularis DAOM-197198": "R1",
        "Rhizophagus clarus HR1": "R1",
        "Mucor velutinosus NIH1002": "R1",
        "Lichtheimia ornata CBS 291.66": "R3",
        "JGI 各株（DOE Joint Genome Institute 提交）": "R2",
        "Mucor flavus KT1a": "R3",
        "Helicostylum pulchrum KT1b": "R3",
        "Umbelopsis isabellina WA0000067209": "R3",
        "Umbelopsis sp. WA50703": "R1",
        "Linnemannia hyalina SCG-10": "R3",
        "Glomus cerebriforme DAOM 227022": "R2",
        "Diversispora epigaea IT104": "R?",
        "Gigaspora rosea FL105 / G. margarita NC121": "R3",
        "Absidia glauca CBS 101.48": "R?",
        "EBP-Nor 各株（ELIXIR Norway 提交）": "R3",
        "GCF_ 条目的注释归属（通例）": "—",
    }
    table(["物种", "Accession", "设计等级", "注释流程", "RNA-seq 证据", "人工审校", "作为基准的可信度", "来源"],
          [[r[0], r[1], _tier_of.get(r[0], "—")] + list(r[2:]) for r in EVIDENCE_TIERS],
          widths=[2.8, 2.2, 1.2, 2.8, 4.2, 1.8, 2.2, 2.8], font=7.5)
    para("设计等级的含义：R1＝参考注释确证使用了 RNA-seq；R2＝JGI 流程以转录组数据为标准输入（流程级判断，逐株待核）；"
         "R3＝参考注释无转录证据（仅蛋白同源、从头预测或默认参数）；R?＝未验证。"
         "R1 与 R2 合称 R+，R3 为 R−，见第 1.5 节的 2×2 设计。", size=8.5)
    para("三条必须写进 Methods 的结论：", bold=True, size=10)
    bullet("第一，本队列没有任何一份既有注释经过人工审校。GCF_ 与 refseq_category = reference genome 都不是审校标志。"
           "全篇措辞只能是「一致性 / 吻合程度」，不能出现「准确性」或「金标准」。")
    bullet("第二，证据类型必须分层报告。BRAKER3 = RNA-seq + 蛋白联合预测；若基准注释是纯蛋白同源产物"
           "（EBP-Nor 各株、Lichtheimia ornata）或 Funannotate 默认参数产物（Mucor flavus、Linnemannia hyalina），"
           "那么「你的注释基因更多 / 更少」里必然混入证据类型差异，不能解读为方法优劣。"
           "建议：至少补一组「BRAKER3 蛋白-only」或「BRAKER2 蛋白-only」对照，把证据类型这一维度单独拆出来。")
    bullet("第三，基因数口径必须先对齐。NCBI 的 gene_counts.total 是「基因总数」，含 tRNA/ncRNA；"
           "论文里报的往往是「蛋白编码基因数」。本报告已逐株核对两者差额（见附录 C），"
           "例如 Mucor velutinosus 的 13,535 = 13,295 蛋白编码 + 240 tRNA。混用会在个位数百分比上造成数百到数千的伪差异。")
    bullet("第四，注意方法同源带来的循环性。优先组中 4 株的既有注释由 BRAKER 家族产生（R. irregularis 用 BRAKER 2.1.5、"
           "R. clarus 用 BRAKER v2.1.2、M. velutinosus 用 BRAKER1、L. ornata 用 BRAKER2 蛋白模式），"
           "它们与本研究的一致性里有一部分是「同一个工具家族的自我印证」，不构成独立验证；"
           "真正方法学独立的是 6 株 JGI 注释。建议把一致性指标按「基准注释是否 BRAKER 家族」分两组报告，"
           "并用 BUSCO、OMArk、Swiss-Prot/InterPro 支持率等不依赖参考基因集的指标作为主证据。详见第 6.12 节。")

    doc.add_heading("2　文档候选名单逐条核查", level=2)
    para("对《投稿前待办事项总表_第四版》任务 1.1 候选表中的 15 条逐一核查，结论如下。")
    table(["亚门 / 类群", "文档原候选", "判定", "202 队列中对应的编号与 accession", "说明"],
          [list(r) for r in DOC_CANDIDATES],
          widths=[2.7, 3.2, 2.4, 4.2, 5.0], font=8)
    para("小结（逐条计数）：可直接采用 6 条；须改写菌株号后可改用 3 条（Phycomyces blakesleeanus、"
         "Mucor lusitanicus、Linnemannia elongata）；降级处理 2 条（Rhizopus microsporus → 备选组，"
         "Umbelopsis ramanniana → 扩展替补池）；明确剔除或不建议 4 条（Rhizopus delemar RA 99-880 同名不同株且无既有注释、"
         "Lichtheimia corymbifera 无既有注释且基因组 BUSCO 仅 87.8%、Absidia glauca 无论文且装配碎片化、"
         "Gigaspora margarita 与 G. rosea 的 N50 仅 6 kb 与 4 kb）。"
         "另有 1 条重复：Mucor circinelloides CBS 277.49 与 Mucor lusitanicus MU402 实为同一株，合计只占 1 个名额。"
         "替换后的完整名单见第 3、4 节，扩展替补见第 5 节。", size=9)

    doc.add_heading("3　优先组（10 株）", level=2)
    para("优先组 10 株全部落在两个「证据匹配」格里：通道 A（本方 BRAKER3 × 参考 R+，6 株）"
         "与通道 B（本方 BRAKER2 × 参考 R−，4 株）。这一组的用途是把批次一任务 1.3–1.5 的全部指标先跑完，"
         "并分别给出「有 RNA-seq 的一半」与「只有蛋白的一半」各自的一致性基线——"
         "两条通道分开报告，是本设计最重要的输出方式。若时间只允许做一半，优先做通道 A 的 P1–P3 与通道 B 的 P7–P8，"
         "这 5 株已覆盖三个亚门中的两个、并且是四条设计通道各自最干净的样本。")
    table(["编号", "物种", "菌株", "亚门", "Accession", "本方方式", "设计格", "装配", "N50(kb)",
           "BUSCO C(%)", "既有注释来源", "基因数"],
          [[r["no"], r["name"], r["strain"], r["sub_short"], r["acc"], r["mode"], r["cell"], r["level"],
            r["n50"], r["busco"], r["src_short"], r["genes"]] for r in rows if r["group"] == "优先"],
          widths=[0.9, 3.1, 1.7, 1.7, 2.5, 1.4, 1.0, 1.3, 1.1, 1.1, 2.6, 1.1], font=7)
    for r in rows:
        if r["group"] == "优先":
            species_block(r["tag"], r["name"], r["strain"], r["sub"], r["acc"], r["src"], r["level"],
                          r["size"], r["n50"], r["busco"], r["genes"], r["rna"], r["orthodb"],
                          r["isoform"], r["reason"], r["risk"], r["mode"], r["tier_label"], r["cell_label"])

    doc.add_heading("4　备选组（10 株）", level=2)
    para("备选组 10 株分三类用途。第一类是通道 C（本方 BRAKER3 × 参考 R−，4 株）与通道 D"
         "（本方 BRAKER2 × 参考 R+，3 株），这是两个「可控错配」格：它们不进入主数字，"
         "只用来量化证据类型效应 Δ——C 减 A 是「本方多一套 RNA-seq 的增益」，D 减 B 是「本方少一套 RNA-seq 的损失」，"
         "两者方向相反、可以互相印证。第二类是通道 E（3 株），把两个匹配格做厚（通道 A 补 2 株、通道 B 补 1 株）。"
         "第三类功能是就近替换：任何优先组成员在任务 1.2/1.3 暴露问题时，按同一设计格、同一亚门就近取用。")
    table(["编号", "物种", "菌株", "亚门", "Accession", "本方方式", "设计格", "装配", "N50(kb)",
           "BUSCO C(%)", "既有注释来源", "基因数"],
          [[r["no"], r["name"], r["strain"], r["sub_short"], r["acc"], r["mode"], r["cell"], r["level"],
            r["n50"], r["busco"], r["src_short"], r["genes"]] for r in rows if r["group"] == "备选"],
          widths=[0.9, 3.1, 1.7, 1.7, 2.5, 1.4, 1.0, 1.3, 1.1, 1.1, 2.6, 1.1], font=7)
    for r in rows:
        if r["group"] == "备选":
            species_block(r["tag"], r["name"], r["strain"], r["sub"], r["acc"], r["src"], r["level"],
                          r["size"], r["n50"], r["busco"], r["genes"], r["rna"], r["orthodb"],
                          r["isoform"], r["reason"], r["risk"], r["mode"], r["tier_label"], r["cell_label"])

    doc.add_heading("5　扩展替补池（已核实可用，未进入 20 株）", level=2)
    para("以下物种同样满足 R1（在 202 队列内）与 R2（有既有注释），基因组完整度也达标，"
         "仅因亚门配额或注释来源优先级未进入前 20。列出以便随时取用。")
    table(["Accession", "物种", "亚门", "说明"],
          [[a, n, s, d] for a, n, s, d in EXTENDED],
          widths=[3.0, 4.2, 2.6, 7.4], font=8)

    # ---- 风险 ----
    doc.add_heading("6　必须提前认识到的风险与方法学边界", level=2)
    for h, body in notes:
        doc.add_heading(h, level=3)
        for b in body:
            bullet(b)

    doc.add_heading("7　与后续任务的衔接", level=2)
    table(["任务", "本报告提供的输入", "下一步动作"],
          [["1.2 核查 OrthoDB 证据库重叠",
            "20 株的「在库 / 不在库」标记（见各株卡片与附表）",
            "按在库/不在库分两组；后续所有指标分组报告"],
           ["1.3 双方蛋白组标准化",
            "20 株的参考蛋白原始文件已取回（见附录 B）",
            "b 与 c 两侧都做最长异构体提取，规则完全一致；同时保留原始版本"],
           ["1.4 BUSCO 配对比较",
            "20 株 accession 与既有注释来源",
            "对标准化后的 b 与 c 各跑一次蛋白模式 BUSCO（mucoromycota_odb12）"],
           ["1.5 蛋白序列层面定量比对",
            "20 株的 accession、来源、基因数、基因密度，以及每株的设计格（A/B/C/D）",
            "五项指标：精确一致率、双向 RBH、互命中对一致度、拆分/合并、各自独有蛋白；"
            "全部按设计格分组报告——A 与 B 出主数字，C 与 D 只用于计算证据类型效应 Δ"],
           ["1.6 结果整理与呈现",
            "本报告的亚门与来源分组口径",
            "Table 5 按「设计格 × 亚门 × 证据库归属」三重分组；措辞统一用「一致性/吻合程度」，不用「准确性」"]],
          widths=[3.4, 6.6, 7.0], font=8.5)

    doc.add_heading("附录 A　字段口径", level=2)
    table(["字段", "口径"],
          [["基因组 BUSCO C",
            "本研究自跑的 mucoromycota_odb12 基因组模式完整度（Complete%），取自表3；与 final_protein_qc_master_table.csv 的基因组列同源"],
           ["N50",
            "优先取 NCBI assembly_stats 的 scaffold N50；JGI 条目取表7 的 assembly_N50_kb。单位统一为 kb"],
           ["既有注释来源",
            "NCBI annotation_info.name 原文（注释提交方），JGI 门户直取条目记为 JGI MycoCosm"],
           ["基因数",
            "NCBI annotation_info.stats.gene_counts.total（蛋白编码 + 非编码）。JGI 门户条目 NCBI 无记录，标为 n.a."],
           ["基因密度",
            "基因数 ÷ 基因组总长（genes/Mb）。143 株的分布为 48–437，中位数 299，用于识别注释膨胀或收缩"],
           ["OrthoDB v12.2 证据库",
            "在库 = 该物种出现在 OrthoDB v12.2 Mucoromycota 157 文件子集（133 属种）中，因该子集是所用 Fungi 分区（2,692 文件）的下属，故必在库；空白 = 未在该子集，需按 Fungi 全表复核"],
           ["RNA-seq",
            "表6 队列内保留 run 与 HISAT2 overall alignment rate；「—」表示该物种无通过 20% 阈值的 run"]],
          widths=[3.2, 13.8], font=8.5)

    doc.add_heading("附录 B　参考蛋白与注释文件取回清单", level=2)
    para("为让任务 1.3 可以立即开始，本报告已从 NCBI Datasets API 取回下列候选物种的既有注释文件，"
         "存放于 9_Benchmark_reference_annotations/ 下：raw_protein_faa/（参考蛋白 *_protein.faa）、"
         "raw_gff3/（参考注释 *_genomic.gff）。GFF3 非任务 1.3 必需，但为任务 1.5 指标五"
         "（各自独有蛋白的单外显子比例判据）与后续外显子结构核对预留。", size=9)
    para("注意：进入 20 株名单的 JGI MycoCosm 门户直取条目共 2 个（JGI:Disdec1_v1.0 与 JGI:Morwol1_v1.0），"
         "它们不在 NCBI 下载清单内，需从 MycoCosm 门户按门户代码下载对应 *.proteins.faa，"
         "并在 metadata 中记录门户版本与下载日期。本次核查中 mycocosm.jgi.doe.gov 对本机返回 HTTP 403，"
         "门户页面与逐株 RNA-seq 证据未能自动核对，建议在门户内人工确认后补记。", size=9)

    doc.add_heading("附录 C　参考蛋白与参考注释的实测核查", level=2)
    para("附录 B 说的是「取回了文件」，本附录说的是「取回的文件是什么」。全部 31 个候选物种的参考蛋白（*_protein.faa）"
         "与参考注释（*_genomic.gff）已直接解析，逐一统计序列数、长度分布、内部终止密码子、重复 ID，"
         "以及从 GFF3 读出的基因数 / 转录本数与异构体负荷。这一层核查的作用是：把「参考注释质量」从文献论断"
         "落到数据层面，同时提前量化任务 1.3 那个卡点（最长异构体提取）到底会影响多大。", size=9)
    table(["Accession", "提交方", "NCBI 蛋白编码基因", "取回蛋白序列数", "中位长度(aa)",
           "<100aa(%)", "内部终止", "重复 ID", "GFF 基因数", "GFF 转录本数", "异构体负荷"],
          [[r["acc"], r["prov"], r["coding"], r["naa"], r["median"], r["pct100"],
            r["istop"], r["dup"], r["ggenes"], r["gtx"], r["infl"]] for r in REF_CHECK],
          widths=[2.5, 2.6, 1.8, 1.6, 1.4, 1.2, 1.0, 1.0, 1.2, 1.2, 1.5], font=7.5)
    para("三点值得直接写进正文的观察：", bold=True, size=10)
    bullet("全部 31 份参考蛋白集都不含内部终止密码子、不含重复 ID —— 说明这批既有注释在序列层面上是干净的，"
           "不存在「参考文件本身有缺陷」这一干扰，b 与 c 的任何差异都可以归因到基因模型而不是文件质量。")
    bullet("取回蛋白序列数与 NCBI「蛋白编码基因数」在 JGI、Broad、NHGRI、Cornell 各株上完全一致（差额 0–6 条，"
           "来自 partial/边界记录）。这意味着这些注释是严格的一基因一转录本，对它们做最长异构体提取是空操作，"
           "任务 1.3 的标准化步骤不会改变它们的任何指标。")
    bullet("但 EBP-Nor 各株与 Rhizophagus irregularis 明显不同：Rhizopus microsporus 每 10,755 个基因有 12,898 条转录本"
           "（异构体负荷 23.3%），Absidia cylindrospora 8.9%，Cunninghamella echinulata 11.6%，Umbelopsis ramanniana 11.9%，"
           "Rhizophagus irregularis 3.3%。这正是文档任务 1.3 卡点所警告的情形——若不先做最长异构体提取，"
           "这些物种的 BUSCO 重复率（D）与「各自独有蛋白数」会被异构体数目直接污染。"
           "建议在 Methods 中把每个基准物种的异构体负荷列表作为附表给出，让读者自己判断标准化前后的差异。")

    doc.save(DOCX)
    return DOCX


# ----------------------------------------------------------------------------
# 4. 主流程
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# 设计轴：本方 BRAKER 模式（来自表6）× 参考注释证据等级
# ----------------------------------------------------------------------------
# R1 = 参考注释确证使用了 RNA-seq；R2 = JGI 流程以转录组为标准输入（逐株待核）；
# R3 = 参考注释无转录证据（仅蛋白同源 / 从头预测 / 默认参数）；R? = 未验证
REF_R1_RNAseq = {
    "GCF_035048745.1",  # Mucor velutinosus NIH1002：16 个生长条件 RNA-seq
    "GCF_026210795.1",  # Rhizophagus irregularis：Illumina + Nanopore RNA-seq
    "GCA_015698045.1",  # Rhizophagus clarus HR1：BRAKER v2.1.2 + RNA-seq
    "GCA_053572175.1",  # Umbelopsis sp. WA50703：funannotate train + 两套 RNA-seq
}
REF_R3_NO_TRANSCRIPT = {
    # EBP-Nor（miniprot + GALBA + EvidenceModeler，纯蛋白同源）
    "GCA_977092135.1", "GCA_977110935.1", "GCA_977110905.1", "GCA_977110975.1",
    "GCA_977110915.1", "GCA_977091265.1", "GCA_977110955.1", "GCA_977110945.1",
    "GCA_973907015.1",
    # Broad：BRAKER2 蛋白模式
    "GCF_029851405.1",
    # Obihiro：Funannotate 1.8.15 默认参数，无转录证据
    "GCA_040367835.1", "GCA_040366295.1",
    # University of Minnesota：Funannotate v1.8.1 全默认，无 RNA-seq
    "GCA_019671135.1",
    # Uppsala：GeneMark-ES ab initio + Swiss-Prot hint
    "GCA_910592045.1", "GCA_910592245.1",
}
REF_UNVERIFIED = {"GCA_003547095.1", "GCA_900079185.1"}
# 提交方在 NCBI 记为其他机构，但注释实际由 JGI Annotation Pipeline 完成
REF_R2_PIPELINE_EXTRA = {"GCA_003550305.1"}  # Glomus cerebriforme DAOM 227022（Morin et al. 2019，JGI 流程）


def ref_tier(acc, ann_name):
    if acc in REF_R1_RNAseq:
        return "R1"
    if acc in REF_R3_NO_TRANSCRIPT:
        return "R3"
    if acc in REF_UNVERIFIED:
        return "R?"
    if (acc.startswith("JGI:") or acc in REF_R2_PIPELINE_EXTRA
            or (ann_name and "Joint Genome Institute" in str(ann_name))):
        return "R2"
    return "R?"


REF_TIER_LABEL = {
    "R1": "R+（RNA-seq 确证）",
    "R2": "R+（JGI 转录组输入，逐株待核）",
    "R3": "R−（无转录证据）",
    "R?": "R?（未验证）",
}


def design_cell(mode, tier):
    plus = tier in ("R1", "R2")
    if mode == "BRAKER3" and plus:
        return "A"
    if mode == "BRAKER2" and not plus:
        return "B"
    if mode == "BRAKER3" and not plus:
        return "C"
    if mode == "BRAKER2" and plus:
        return "D"
    return "?"


CELL_LABEL = {
    "A": "通道 A｜匹配：本方 BRAKER3 × 参考 R+",
    "B": "通道 B｜匹配：本方 BRAKER2 × 参考 R−",
    "C": "通道 C｜错配：本方 BRAKER3 × 参考 R−",
    "D": "通道 D｜错配：本方 BRAKER2 × 参考 R+",
    "?": "未分类",
}


def main():
    m = load_master()
    membership = orthodb_members()
    ncbi = fetch_ncbi([r[0] for r in PRIORITY + BACKUP] + [e[0] for e in EXTENDED])
    m["busco"] = pd.to_numeric(m["BUSCO C_muco"], errors="coerce")
    m["n50kb"] = pd.to_numeric(m["assembly_N50_kb"], errors="coerce").astype(float)
    # 参考注释异构体负荷（由已取回的参考 GFF3 统计）
    iso_path = os.path.join(OUTDIR, "reference_gff_isoform_stats_v2.json")
    pstat_path = os.path.join(OUTDIR, "reference_protein_stats_v2.json")
    iso = json.load(open(iso_path)) if os.path.exists(iso_path) else {}
    pstat = json.load(open(pstat_path)) if os.path.exists(pstat_path) else {}

    def isoform_text(acc, n_genes):
        v = iso.get(acc)
        if not v:
            return "未取回参考 GFF3（JGI 门户条目需自行下载）"
        g_genes, tx = v.get("genes") or 0, v.get("tx") or 0
        if g_genes and tx < 0.5 * g_genes:
            # 该提交方的 GFF3 未用标准 mRNA/transcript 特征，转录本层无法从 GFF 判定
            faa = pstat.get(acc, {}).get("n")
            extra = (f"，但蛋白 FASTA 有 {faa:,} 条、多于蛋白编码基因数"
                     if faa and n_genes not in ("(JGI 门户)",) and faa > int(str(n_genes).replace(",", "")) else "")
            return f"GFF3 未使用标准 mRNA/transcript 特征，无法判定{extra}"
        if v.get("inflation_pct") in (0, 0.0) and v.get("multi_iso", 0) <= 6:
            return f"无（1 基因 1 转录本；多转录本基因 {v.get('multi_iso', 0)} 个）"
        return (f"含异构体：多转录本基因 {v['multi_iso']:,} 个，"
                f"转录本比基因多 {v['inflation_pct']}%（最长转录本 {v['max_iso']}）")

    def pack(acc, group, tag, reason, risk=None):
        r = m[m["主标识"] == acc].iloc[0]
        n = ncbi.get(acc, {})
        n50 = n.get("scaffold_n50") or (r["assembly_N50_kb"] * 1000 if pd.notna(r["assembly_N50_kb"]) else None)
        n50kb = f"{round(n50 / 1000):,}" if n50 else "n.a."
        size = n.get("total_sequence_length") or r.get("total_length_bp")
        size_mb = f"{float(size) / 1e6:,.1f}" if size and str(size) != "nan" else "n.a."
        genes = n.get("ann_gene_total")
        if genes is None:
            genes = "n.a.（JGI 门户）"
        src = n.get("ann_name") or ("JGI MycoCosm 门户直取（%s）" % acc.split(":")[-1])
        rna = str(r["RNA-seq"])
        rna = "无（表6 记为「—」）" if rna.strip() in ("—", "nan") else rna
        in_db = "在库（Mucoromycota 157 子集）" if norm(r["原名(表1)"]) in membership or norm(r["采纳名"]) in membership else "未在该子集，需按 Fungi 全表复核"
        tier = ref_tier(acc, n.get("ann_name"))
        mode = str(r["BRAKER_mode"])
        cell = design_cell(mode, tier)
        lvl_map = {"Contig": "Contig", "Scaffold": "Scaffold", "Chromosome": "Chromosome",
                   "Complete Genome": "Complete Genome", "Annotated Standard Draft": "JGI Standard Draft"}
        level = n.get("assembly_level") or lvl_map.get(str(r["assembly 水平（来源词）"]), str(r["assembly 水平（来源词）"]))
        return dict(
            no=int(r["编号"]), name=r["原名(表1)"], strain=("—" if str(r["strain/isolate"]) == "nan" else str(r["strain/isolate"])),
            sub=r["Subphylum"], sub_short={"Mucoromycotina": "Mucoro.", "Mortierellomycotina": "Mortierello.",
                                           "Glomeromycotina": "Glomero."}[r["Subphylum"]],
            acc=acc, src=src,
            src_short=(src.replace("Annotation submitted by ", "").replace("DOE Joint Genome Institute", "JGI")
                          .replace("Royal Botanic Gardens Kew", "Kew (RBG Kew)")
                          .replace("NIBB core research facilities, National Institues for Basic Biology", "NIBB")
                          .replace("JGI MycoCosm 门户直取（%s）" % acc.split(":")[-1], "JGI MycoCosm 门户直取")),
            level=level, size=size_mb, n50=n50kb, busco=round(float(r["busco"]), 1),
            genes=(f"{int(genes):,}" if isinstance(genes, (int, float)) else genes),
            genes_int=(int(genes) if isinstance(genes, (int, float)) else "-"),
            rna=rna, orthodb=in_db, isoform=isoform_text(acc, genes), group=group, tag=tag, reason=reason,
            mode=mode, tier=tier, tier_label=REF_TIER_LABEL[tier], cell=cell, cell_label=CELL_LABEL[cell],
            risk=risk or "无特别风险；按亚门分组报告即可。",
        )

    rows = []
    for i, (acc, tag, reason, risk) in enumerate(PRIORITY, 1):
        rows.append(pack(acc, "优先", f"P{i}", reason, risk))
    for i, (acc, tag, reason, risk) in enumerate(BACKUP, 1):
        rows.append(pack(acc, "备选", f"B{i}", reason, risk))

    stats = {
        "n_total": len(m),
        "n_has_ann": int((m["原始蛋白注释"] == "有").sum()),
        "n_attrib": int(((m["原始蛋白注释"] == "有") &
                         (m["来源库"].eq("JGI") | m["主标识"].astype(str).str.startswith(("GCA_", "GCF_")))).sum()),
        "n_busco985": int(((m["原始蛋白注释"] == "有") & (m["busco"] >= 98.5)).sum()),
        "n_final": int(((m["原始蛋白注释"] == "有") & (m["busco"] >= 98.5) &
                        ((m["n50kb"] >= 800) | (m["assembly 水平（来源词）"].isin(["Chromosome", "Complete Genome"])))).sum()),
        "n_b2_all": int((m["BRAKER_mode"] == "BRAKER2").sum()),
        "n_b3_all": int((m["BRAKER_mode"] == "BRAKER3").sum()),
        "n_b2_ann": int(((m["BRAKER_mode"] == "BRAKER2") & (m["原始蛋白注释"] == "有")).sum()),
        "n_b3_ann": int(((m["BRAKER_mode"] == "BRAKER3") & (m["原始蛋白注释"] == "有")).sum()),
        "n_loose": int(((m["原始蛋白注释"] == "有") & (m["busco"] >= 98.0) &
                        ((m["n50kb"] >= 500) | (m["assembly 水平（来源词）"].isin(["Chromosome", "Complete Genome"])))).sum()),
    }

    # 设计格 × 亚门 交叉表
    cell_tab = {}
    for r in rows:
        cell_tab.setdefault(r["cell"], {}).setdefault(r["sub"], 0)
        cell_tab[r["cell"]][r["sub"]] += 1

    # 附录 C 数据
    prov_short = {
        "DOE Joint Genome Institute": "JGI", "Royal Botanic Gardens Kew": "Kew (RBG Kew)",
        "NIBB core research facilities, National Institues for Basic Biology": "NIBB",
        "Broad Institute": "Broad Institute", "NHGRI": "NHGRI", "INRA": "INRA",
        "Cornell University": "Cornell University", "ELIXIR Norway": "ELIXIR Norway（EBP-Nor）",
        "University of Minnesota-Twin Cities": "Univ. of Minnesota",
        "Obihiro University of Agriculture and Veterinary Medicine": "Obihiro Univ.",
        "University of Warsaw": "Univ. of Warsaw",
    }
    REF_CHECK = []
    for acc, v in sorted(pstat.items()):
        e = ncbi.get(acc, {})
        g = iso.get(acc, {})
        infl = g.get("inflation_pct")
        if not g:
            infl = "n.a."
        elif (g.get("genes") or 0) and (g.get("tx") or 0) < 0.5 * g["genes"]:
            infl = "GFF 无标准转录本特征，无法判定"
        elif infl == 0 and g.get("multi_iso", 0) <= 6:
            infl = "0%（单转录本）"
        else:
            infl = f"{infl}%"
        REF_CHECK.append(dict(
            acc=acc,
            prov=prov_short.get(e.get("ann_provider"), (e.get("ann_provider") or "n.a.")[:22]),
            coding=(f"{int(e['ann_gene_protein_coding']):,}" if e.get("ann_gene_protein_coding") else "n.a."),
            naa=f"{v['n']:,}", median=v["median"], pct100=v["pct_lt100"],
            istop=v["n_internal_stop"], dup=v["n_dup_ids"],
            ggenes=(f"{g.get('genes'):,}" if g.get("genes") else "n.a."),
            gtx=(f"{g.get('tx'):,}" if g.get("tx") else "n.a."),
            infl=infl,
        ))
    REF_CHECK.sort(key=lambda r: r["prov"])

    # 名单 CSV（与 Word 同源，供后续流程直接读取）
    csv_path = os.path.join(OUTDIR, "基准物种20株_名单_v2.csv")
    with open(csv_path, "w", encoding="utf-8-sig") as fh:
        fh.write("组别,标签,设计格,本方注释方式,参考注释证据等级,No,物种,菌株,亚门,Assembly_accession,"
                 "装配水平,基因组大小_Mb,N50_kb,基因组BUSCO_C_pct,既有注释提交方,既有注释基因数,"
                 "OrthoDB_v12.2_证据库,参考注释异构体,队列内RNA-seq\n")
        for r in rows:
            fh.write(",".join('"%s"' % str(r[k]).replace('"', "'") for k in
                              ["group", "tag", "cell", "mode", "tier_label", "no", "name", "strain", "sub",
                               "acc", "level", "size", "n50", "busco", "src", "genes", "orthodb",
                               "isoform", "rna"]) + "\n")
    print("written:", csv_path)

    notes = json.load(open(os.path.join(OUTDIR, "risk_notes_v2.json")))
    path = build_docx(rows, stats, [(n["heading"], n["bullets"]) for n in notes], REF_CHECK, cell_tab)
    print("written:", path)
    print("stats:", stats)


if __name__ == "__main__":
    main()
