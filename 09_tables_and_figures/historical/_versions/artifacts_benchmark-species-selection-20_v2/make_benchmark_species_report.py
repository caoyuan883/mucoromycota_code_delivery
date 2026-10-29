#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成《基准物种遴选报告》Word 文档（批次一 · 任务 1.1）。

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
CACHE = os.path.join(OUTDIR, "benchmark_species_evidence.json")
DOCX = os.path.join(OUTDIR, "基准物种遴选报告_20株_v1.docx")

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
PRIORITY = [
    ("GCF_026210795.1", "Glomeromycotina 旗舰；全队列唯一染色体级 AMF",
     "三个亚门里 Glomeromycotina 装配质量最差、最易被质疑，本株是全队列 202 株中唯一达到染色体级（Chromosome）的 AMF 装配：32 条染色体、42 contig → 32 scaffold、contig N50 3.90 Mb、scaffold N50 5.09 Mb、仅 10 个 gap，其中 23 条染色体完全无 gap、17 条达到端粒到端粒；mucoromycota_odb12 基因组模式 BUSCO 完整度 C = 99.8%。既有注释由 Royal Botanic Gardens Kew 提交，30,147 个蛋白编码基因、31,153 条转录本，且注释有短读 + 长读双 RNA-seq 证据（Illumina RNA-seq → STAR → BRAKER 2.1.5；Nanopore PCR-cDNA RNA-seq → funannotate 1.8.7 update + PASA 精修 UTR 与 poly(A)），并经人工剔除含转座子结构域的基因。原始论文为 Manley et al. 2023, G3 13(6):jkad077（DOI 10.1093/g3journal/jkad077）。DAOM 197198 是 AMF 中注释版本最多、被引用最广的菌株，以它为基准，可以把「我的注释是否达到 AMF 领域公认水平」回答得最干净；队列内亦保留了两条高比对率 RNA-seq（SRR3285886 98.96%、SRR3285885 98.93%）。",
     "三条必须写进 Methods 的限制。第一，注释流程本身是 BRAKER 2.1.5——与本研究 BRAKER2/3 属同一工具家族，因此「一致」在多大程度上是方法同源造成的，必须显式讨论，不能用它单独证明方法正确。第二，基因数存在三套不一致的版本：论文报 30,209，GenBank 提交版（GCA_026210795.1）为 30,190 个蛋白编码基因，RefSeq 版（GCF_026210795.1）为 30,147 个——本报告锁定的对象是 GCF_026210795.1，引用时必须写明取的是哪一套，不能把论文数字套到提交版上。第三，GCF 记录挂在独立的 RefSeq BioProject PRJNA1128271 下，RefSeq 侧只是重新发布提交方注释（其 GenBank 孪生记录为 GCA_026210795.1），不可写成「NCBI 注释」。此外基因组 146.8 Mb，基因密度仅 205 genes/Mb，明显低于 Mucorales，比较时应按亚门分别报告。"),    ("GCA_015698045.1", "Glomeromycotina 第二代表；与上一株构成属内对照",
     "基因组 BUSCO C = 99.8%，contig N50 = 1.31 Mb（360 个 contig，L50 = 40），是 Glomeromycotina 内除 R. irregularis 之外装配连续性最好的一株，基因组 146.8 Mb 与前者体量相当，构成天然的属内对照。既有注释由日本基础生物学研究所（NIBB）提交（2020-11-13），30,348 个基因（30,194 蛋白编码 + 154 nc），蛋白模式为 BRAKER v2.1.2 + AUGUSTUS 3.3.2 + GeneMark-ES 4.32，作者自述为 RNA-seq-based，即证据类型与 BRAKER3 部分对齐。队列内有两条 RNA-seq（DRR041161 88.24%、DRR041163 85.86%）。BioProject PRJDB8680 的标题即 Improved genome sequencing of arbuscular mycorrhizal fungi, Rhizophagus clarus HR1。",
     "必须注明的三个限制。第一，该基因组的唯一文献是 2026-08-13 发布于 Research Square、状态为 Under Revision 的预印本（DOI 10.21203/rs.3.rs-10340957/v1），未经同行评审；组装自 2020-11 释放后长期没有正式发表。作者在预印本中明确表示该组装未达染色体级连续性。第二，基因数在预印本（31,223 蛋白编码 / 32,541 蛋白模型）与 NCBI（30,348）之间不一致，本报告锁定 NCBI 提交版。第三，该株的参考 GFF3 未使用标准 mRNA/transcript 特征（只有 gene 与 CDS，transcript 特征仅 22 条），无法从 GFF 读出异构体结构，但参考蛋白 FASTA 有 32,541 条，明显多于 30,194 个蛋白编码基因，说明蛋白集含异构体——任务 1.3 的最长异构体提取对这一株必须以蛋白 FASTA 为准。另需区分：2018 年 BMC Genomics 那篇（DOI 10.1186/s12864-018-4853-0）用的是另一套组装（BEXD01，27,753 基因），引用时不可混用。"),    ("GCF_035048745.1", "全部候选中唯一有 RNA-seq 支撑、且证据类型与 BRAKER3 最可比的基准",
     "这是本报告最推荐的一株，理由有三层。其一，证据类型对齐：参考注释使用 16 个生长条件（RPMI ± 人 PMN、YG/MEA 不同温度、铁限制、H2O2、SBA 等；SAMN41019669–84）的 RNA-seq，经 STAR v2.5.3a → StringTie v1.3.3，再与 TransDecoder v5.0.2、BRAKER1 v1.9、CodingQuarry v2.0 的结果一并交给 EVidenceModeler 整合，最后由 funannotate 收尾——这与 BRAKER3「RNA-seq + 蛋白联合」的证据结构最接近，比较结果不会混入'证据类型不同'这一混杂因素。其二，方法学可核实：有专门的论文（Conlan et al. 2024, Microbiol Resour Announc 13(7):e00062-24, DOI 10.1128/mra.00062-24），流程每一步都有版本号。其三，装配好：42 个 contig、contig N50 2.65 Mb、Merqury QV 46.6 / 完整度 99.7%、17 个 contig 带 TTAGGGn 端粒（提示 9 条染色体）。该株为 RefSeq 参考基因组，基因组 BUSCO C = 99.7%，队列内有两条 RNA-seq（90.19%、89.73%）。",
     "基因数口径要当心：NCBI 的 13,535 是「基因总数」，其中蛋白编码 13,295 + tRNA 240；论文报的 13,295 是蛋白编码数，两者完全对账。另有一处小矛盾：NCBI 装配记录的注释说明栏写 Funannotate v1.8.14，而论文 Methods 写 v1.8.13（本报告未能取得该说明栏原文，标注为待核）。表3 采纳名为 Mucor circinelloides、原名为 Mucor velutinosus，是 202 队列中 5 组同物异名之一；写作时必须同时给出原名与采纳名，并说明与 No.120（Mucor lusitanicus）不是同一株。该注释同样无人工审校。另需注意方法同源：其流程包含 BRAKER1 v1.9，属 BRAKER 工具家族，因此本株也带一定循环性（见第 6.10 节）；不过它与 BRAKER3 的差异主要集中在证据组合方式（BRAKER1 只用 RNA-seq，本研究的 BRAKER3 为 RNA-seq + 蛋白），这一点恰恰可以在讨论中正面利用。"),
    ("GCA_040256725.1", "Mucoromycotina 模式菌株；全队列基因组完整度最高一档",
     "mucoromycota_odb12 基因组模式 BUSCO C = 99.9%，与全队列最高值并列；N50 = 4.61 Mb，为 Mucor 属内最好。既有注释为 JGI 提交（Mucci3，2024-06-21），11,759 个基因，密度 321 genes/Mb，落在 Mucor 属常规区间内，无膨胀迹象。队列内两条 RNA-seq 比对率均在 99% 以上（SRR14527209 99.20%、SRR14527207 99.19%），是全部候选物种中转录组证据最强的一株。该株即文档候选表中的'Mucor circinelloides CBS 277.49'——现用名为 Mucor lusitanicus，同一菌株号，命名变更须在表1/表3口径下书写。",
     "文档原候选写作'MU402'，但 202 队列与 NCBI 记录均为 CBS 277.49；MU402 并非本队列菌株，报告中必须改写为 CBS 277.49，否则菌株身份对不上。"),
    ("GCA_040209145.1", "经典模式生物；对应文档候选（菌株号须改写）",
     "基因组 BUSCO C = 99.6%，N50 = 1.89 Mb。既有注释为 JGI 提交（PhyblU21_2，2024-06-17），13,146 个基因，密度 232 genes/Mb。队列内两条 RNA-seq 比对率 97.85%/97.56%（SRR10402501、SRR9002717）。Phycomyces 是光生物学与类胡萝卜素研究的经典模式属，其基因模型在文献中可交叉比对，出现异常时容易定位原因。",
     "文档原候选写作'NRRL1555'，202 队列中实际为 UBC21（GCA_040209145.1），两者为不同菌株，必须改写。该基因组偏大（56.7 Mb）而基因密度偏低（232 genes/Mb），属高重复含量装配，解释指标时应说明。"),
    ("GCF_025201355.1", "RefSeq 参考基因组 + JGI 注释 + 高装配连续性，三重占优",
     "同时满足文档判据的第 1 条与第 3 条：accession 为 GCF_ 开头的 RefSeq 参考基因组（refseq_category = reference genome），且 N50 = 3.62 Mb 为全队列第 4 高、Mucorales 内前 3；基因组 BUSCO C = 99.6%。既有注释由 JGI 完成（Halrad1，2022-11-01），10,305 个基因，密度 393 genes/Mb；队列内 RNA-seq SRR7476973 比对率 98.72%。基因组仅 26.2 Mb 却达到 3.62 Mb N50，是'小而完整'的典型，适合作为装配质量与注释质量解耦的样本。",
     "基因密度 393 genes/Mb 位于 143 株分布的上四分位以上，需在结果中与队列分布并列展示，避免被误读为注释膨胀。"),
    ("GCF_029851405.1", "Lichtheimia 属唯一可用者；定位为「方法学对照」而非外部金标准",
     "文档候选表中的 Lichtheimia corymbifera（No.75，GCA_037042095.1）在本队列中'原始蛋白注释 = 无'，且基因组 BUSCO 仅 87.8%，两条判据均不满足，必须替换。同属的 Lichtheimia ornata（No.77）为 GCF_ 开头的 RefSeq 参考基因组，BUSCO C = 99.8%，既有注释由 Broad Institute 提交（2023-07-27），13,169 个基因，是 Lichtheimia 属（医学重要类群，接合菌病主要病原之一）在本队列中唯一既有注释可用者。该注释的原始论文为 Shea et al. 2023, Microbiol Resour Announc 12(7):e00338-23（DOI 10.1128/mra.00338-23）。需要特别指出：论文 Methods 写明该注释由 BRAKER2 蛋白模式产生，证据为来自 L. corymbifera / L. ramosa / L. hyalospora 的 35,966 条蛋白，未使用 RNA-seq。因此它与本研究 BRAKER2 那一层的注释属于同一工具家族，严格说是「方法学对照」而不是独立的外部金标准；这一点必须写进 Methods，不能含糊。",
     "该株在队列内无保留 RNA-seq run（表6 记为'—'），参考注释本身亦为 BRAKER2 蛋白模式、无转录证据 → 它是一套「BRAKER2 蛋白-only 基线」，不是金标准。装配为 Scaffold 级、N50 = 386 kb、contig N50 仅 151 kb（ONT + Illumina + TELL-Seq，Supernova 组装），碎片化明显。以蛋白层面一致性比较为目的可以接受，但不得用于外显子边界或装配质量论证，也不得表述为「与人工审校注释比较」。"),
    ("JGI:Disdec1_v1.0", "全队列装配连续性最好的一株；Mortierellomycotina 代表",
     "N50 = 6.89 Mb，为 202 株中最高值；基因组 BUSCO C = 99.9%，与全队列最高值并列。既有注释由 JGI 完成，Mortierellomycotina 内稀缺的高质量装配。以它作为 Mortierellomycotina 的第一个基准，可以在'基因组质量最好'的条件下先确立注释质量的基线，再去看该亚门其余较破碎的装配。",
     "属 JGI MycoCosm 门户直取（非 NCBI accession），比对时须同步记录 JGI 门户代码与版本号 Disdec1_v1.0，仓储 metadata 中要写清下载日期与 URL。队列内有 RNA-seq run SRR13647699（97.75%）可用于 BRAKER3 侧；参考注释本身的转录组证据仍需回到 JGI 门户核实（MycoCosm 对外访问返回 403，建议在门户内人工核对）。"),
    ("GCA_036320915.1", "Mortierellaceae 分类修订的核心材料；对应文档候选（菌株号须改写）",
     "基因组 BUSCO C = 99.8%，N50 = 2.74 Mb，Mortierellomycotina 内仅次于 Disdec1。既有注释为 JGI 提交（MoeNVP64_1，2024-01-31），13,991 个基因，密度 274 genes/Mb。Linnemannia elongata 是 Mortierellaceae 分类修订中的关键物种，注释质量直接影响该科基因家族的解读。队列内两条 RNA-seq（92.77%、90.85%）。与 Disdec1 同科不同属，构成 Mortierellomycotina 内部的属级对照。",
     "文档原候选写作'AG-77'，202 队列实际为 NVP64（GCA_036320915.1），菌株号必须改写。"),
    ("GCF_002105155.1", "Mortierellomycotina 早期分支；发布时间最早、修订轮次最多的 JGI 参考",
     "RefSeq 参考基因组，既有注释由 JGI 完成（Lobtra1，2017-08-29），11,983 个基因，密度 280 genes/Mb，基因组 BUSCO C = 98.6%。该株在 Mortierellomycotina 中处于较早分出的位置，与 Disdec1、Linnemannia elongata 的系统发育距离足够，能够把'注释一致性是否随系统发育距离衰减'这一混杂因素纳入设计。注释自 2017 年发布后经 JGI 多轮更新，是队列中'经过时间检验'的少数注释之一。队列内有 RNA-seq SRR8840873（97.60%）。",
     "N50 = 673 kb，装配连续性中等；GC 与重复含量属 Mortierellomycotina 典型值，比较时按亚门分组。该 GCF 记录在 NCBI 存在配对的 GenBank 孪生记录 GCA_002105155.1（assembly name 均为 Lobtra1），RefSeq 侧只做重新发布，因此同样不可写成「NCBI 注释」。"),
]

BACKUP = [
    ("GCA_002105175.1", "Mucoromycotina",
     "RefSeq 类别为 reference genome 的 GCA 装配，JGI 注释（Absrep1，2017-04-20），15,151 个基因；BUSCO C = 99.4%，N50 = 1.30 Mb；队列内 RNA-seq SRR7974509（96.96%）。JGI 早期参考基因组，注释发布时间早、被复用次数多，且参考注释为单转录本（无异构体）。"),
    ("GCA_002105135.1", "Mucoromycotina",
     "JGI 注释（PRJNA330704），11,247 个基因；BUSCO C = 99.6%，N50 = 2.37 Mb；RNA-seq SRR8843209（98.31%）；参考注释单转录本。与 Absidia repens 同为 JGI 早期 Mucorales 参考批次，便于做批次内一致性检查。"),
    ("GCA_025094155.1", "Mucoromycotina",
     "JGI 注释，12,174 个基因；BUSCO C = 99.7%，N50 = 129 kb；参考注释单转录本。与优先组 P7（Lichtheimia ornata，Broad/BRAKER2）构成同属双基准：一株是 BRAKER2 蛋白-only 产物，一株是 JGI 机构流程产物，可直接回答「同属之内，注释差异来自工具还是来自数据」。注意该株在队列内无保留 RNA-seq run（表6 记为「—」），比较时两层都会是 BRAKER2，需在 Methods 中说明。"),
    ("GCF_025331425.1", "Mucoromycotina",
     "RefSeq 参考基因组，JGI 注释（Radspe1，2022-11-01），11,044 个基因；BUSCO C = 99.2%，N50 = 1.93 Mb；RNA-seq SRR6047893（98.91%）；参考注释单转录本。Radiomyces 属在 Mucorales 内位置较特殊，可补充系统发育采样。"),
    ("GCF_025528875.1", "Mucoromycotina",
     "RefSeq 参考基因组，JGI 注释（Mycafr1，2023-02-06），11,085 个基因；BUSCO C = 99.0%，N50 = 1.96 Mb；RNA-seq SRR6049697（98.68%）；参考注释单转录本。Mycotypha 属代表，与扩展池中的 Mycotypha indica（JGI 门户）构成属内双株，可交叉验证。"),
    ("GCF_025118155.1", "Mucoromycotina",
     "RefSeq 参考基因组，JGI 注释（Cokrec1，2022-11-01），11,066 个基因；BUSCO C = 98.4%，N50 = 811 kb；RNA-seq SRR6056974（91.60%）；参考注释单转录本。基因组仅 28.1 Mb，适合做小基因组对照。"),
    ("GCA_040367835.1", "Mucoromycotina",
     "本组唯一「流程版本 + 论文」双可考的候选：PacBio HiFi（Sequel IIe）+ IPA v1.8.0 组装，69 contig、N50 858 kb；注释用 Funannotate（Galaxy 1.8.15+galaxy1），仅改 --busco_db mucolares 与 --busco_seed_species Rhizopus_oryzae 两项；论文 Hosono et al. 2024, Microbiol Resour Announc 13(9):e00343-24（DOI 10.1128/mra.00343-24）。BUSCO C = 99.8%，11,530 个基因（11,295 蛋白编码 + 235 nc）。无 RNA-seq。"),
    ("JGI:Morwol1_v1.0", "Mortierellomycotina",
     "JGI 门户直取，N50 = 4.65 Mb（全队列第 3 高），BUSCO C = 99.2%。Mortierellomycotina 内装配质量仅次于 Disdec1，可作为该亚门第四株基准。需从 MycoCosm 门户下载参考蛋白并记录门户版本与下载日期。"),
    ("GCA_003547095.1", "Glomeromycotina",
     "Cornell University 提交的 Diversispora epigaea IT104，28,348 个基因，BUSCO C = 98.5%，N50 = 453 kb；参考注释单转录本（1 基因 1 转录本，28,348 : 28,348 完全对齐）。Glomeromycotina 内第三株可选基准，属级多样性与前两株（Rhizophagus）不同；装配较破碎（156.6 Mb），宜作补充而非常规基准。该株在队列内亦无保留 RNA-seq run，本项目对它是 BRAKER2 运行。更关键的限制是：其注释流程与证据基础本次未能核实（原文受出版商访问限制，无可用 OA 版本），属于「来源不明的基因集」，因此只宜作补充。"),
    ("GCA_977110955.1", "Mucoromycotina",
     "文档候选中的 Rhizopus microsporus（临床与工业双重要）。EBP-Nor 装配，BUSCO C = 99.4%，N50 = 4.72 Mb。但参考注释为 EBP-Nor 自研流程（miniprot + GALBA + EvidenceModeler）纯蛋白同源证据的自动产物：无 RNA-seq、无 RefSeq 镜像、无注释方法论文，且异构体负荷高达 23.3%（1,696/10,755 基因有多转录本）。列入备选组是作为「蛋白-only 自动流程基线」的对照样本，不得当作高可信度基准；这一条同时解释了为什么文档把它列在候选表里却不宜进优先组。"),
]

EXTENDED = [
    ("GCA_977092135.1", "Absidia cylindrospora", "Mucoromycotina", "EBP-Nor，BUSCO C=99.4%，N50=6.36 Mb，15,992 基因。装配极好但注释为纯蛋白同源自动流程、异构体负荷 8.9%，定位为「装配好 / 注释弱」的对照。"),
    ("GCA_977110915.1", "Cunninghamella echinulata", "Mucoromycotina", "EBP-Nor，BUSCO C=99.4%，N50=4.26 Mb，11,249 基因，异构体负荷 11.6%；文档候选之一。"),
    ("GCA_977110945.1", "Umbelopsis ramanniana", "Mucoromycotina", "EBP-Nor，BUSCO C=99.4%，N50=1.65 Mb。既有注释仅 9,709 个基因（390 genes/Mb），文档所述「基因数异常」经核实成立；建议作为异常例证用于批次二任务 2.4，不作基准。异构体负荷 11.9%。"),
    ("GCA_977091265.1", "Mortierella alpina", "Mortierellomycotina", "EBP-Nor，BUSCO C=99.3%，N50=3.23 Mb，11,632 基因；工业重要（ARA 生产）。异构体负荷 7.7%。"),
    ("GCA_025201345.1", "Thamnidium elegans", "Mucoromycotina", "JGI，BUSCO C=99.5%，N50=916 kb，13,308 基因；参考注释单转录本。"),
    ("JGI:Abscae1_v1.0", "Absidia caerulea", "Mucoromycotina", "JGI 门户，BUSCO C=99.4%，N50=3.62 Mb。"),
    ("JGI:FenlinCB54_1_v1.0", "Fennellomyces linderi", "Mucoromycotina", "JGI 门户，BUSCO C=99.9%，N50=3.56 Mb。"),
    ("JGI:Benpoi1_v1.0", "Benjaminiella poitrasii", "Mucoromycotina", "JGI 门户，BUSCO C=99.2%，N50=3.50 Mb。"),
    ("JGI:Ellano2_v2.0", "Ellisomyces anomalus", "Mucoromycotina", "JGI 门户，BUSCO C=99.8%，N50=2.29 Mb。"),
    ("JGI:Moramb1_v1.0", "Mortierella ambigua", "Mortierellomycotina", "JGI 门户，BUSCO C=98.6%，N50=3.39 Mb；Mortierellomycotina 第五备选。"),
    ("GCA_040366295.1", "Helicostylum pulchrum", "Mucoromycotina", "PacBio HiFi + IPA，Funannotate 1.8.15，无 RNA-seq；BUSCO C=99.7%，N50=946 kb，12,199 基因；论文同 KT1a。"),
    ("GCA_016758905.1", "Umbelopsis isabellina", "Mucoromycotina", "University of Warsaw；MAKER v2.31.8（AUGUSTUS + GeneMark-ES）用 100,082 条公共 EST + 273,458 条 Mucoromycotina 蛋白训练（是 EST，不是 RNA-seq）；BUSCO C=99.6%，N50=1.34 Mb，9,033 基因；论文 Muszewska et al. 2021, Front Microbiol 12:636986。"),
    ("GCA_019671135.1", "Linnemannia hyalina", "Mortierellomycotina", "University of Minnesota；Funannotate v1.8.1 全默认参数，无 RNA-seq；BUSCO C=99.4%，N50=2.32 Mb，13,010 基因（12,693 蛋白编码 + 317 tRNA）。基因组 51.6 Mb 偏大、基因密度仅 246/Mb（本批最低），用前建议先查重复含量与污染。"),
    ("GCA_053572175.1", "Umbelopsis sp. WA50703", "Mucoromycotina", "University of Warsaw；Funannotate v1.8.15 且用两套 RNA-seq 训练；BUSCO C=99.2%，N50=980 kb。注意同株另有 GCA_964291815.1（8,693 基因），文献里的 MAISEN 人工审校属于后者，不能算在本装配上。"),
    ("GCF_041956525.1", "Rhizomucor pusillus", "Mucoromycotina", "RefSeq + JGI，BUSCO C=98.8%，11,084 基因；临床重要。"),
    ("GCF_025094135.1", "Mucor mucedo", "Mucoromycotina", "RefSeq + JGI，BUSCO C=97.1%，14,226 基因。"),
    ("GCF_025766255.1", "Zychaea mexicana", "Mucoromycotina", "RefSeq + JGI，BUSCO C=98.1%，14,134 基因。"),
    ("GCA_003550305.1", "Glomus cerebriforme", "Glomeromycotina", "INRA；BUSCO C=98.7%，21,336 基因。N50=148 kb，装配破碎，作 AMF 第三选择。"),
    ("GCA_910592045.1", "Gigaspora rosea", "Glomeromycotina", "Uppsala，BUSCO C=98.5%，但 N50 仅 4 kb，装配高度破碎，不建议作基准。"),
    ("GCA_910592245.1", "Gigaspora margarita", "Glomeromycotina", "Uppsala，BUSCO C=98.8%，但 N50 仅 6 kb，装配高度破碎，不建议作基准。"),
    ("GCA_900079185.1", "Absidia glauca", "Mucoromycotina", "不建议采用：未检索到任何发表论文（NCBI assembly→PubMed 链接为空），注释流程无任何文档，Ensembl Fungi 确证为原样导入（genebuild_method = import, 2016-04）；装配仅 Illumina + Velvet，5,338 scaffold、contig N50 仅 11 kb，是候选中最碎片化的。RNA-seq 证据与人工审校情况均无法核实。"),
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

def build_docx(rows, stats, notes, REF_CHECK):
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
                      isoform, reason, risk):
        doc.add_heading(f"{tag}　{name} {strain}", level=3)
        table(
            ["项目", "内容", "项目", "内容"],
            [
                ["亚门", sub, "Assembly accession", acc],
                ["既有注释提交方", src, "装配水平", level],
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
    doc.add_heading("基准物种遴选报告", level=0)
    para("批次一 · 任务 1.1 —— 与高质量参考注释做蛋白层面对照的物种名单", bold=True, size=12, align="c")
    para("Mucoromycota 202 基因组注释项目（Scientific Data Data Descriptor）", size=10, align="c")
    para("版本 v1 · 2026-09-14 · 全部结论以 Assembly accession 锚定，可逐条复核",
         size=9, align="c", space_after=14)

    doc.add_heading("0　一页速览", level=2)
    para(f"从本研究 202 株队列中，按'基因组高质量 + 既有注释蛋白高质量'两条并行的判据，"
         f"筛出 {stats['n_has_ann']} 株具备第三方既有注释的物种，再施加基因组完整度与装配连续性门槛，"
         f"最终确定优先组 10 株、备选组 10 株，共 20 株。全部 20 株均在本研究 202 株队列内，"
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
    bullet("优先组覆盖 3 个亚门：Mucoromycotina 5 株、Mortierellomycotina 3 株、Glomeromycotina 2 株，"
           "与队列 123 : 57 : 22 的结构大致同比例，并对 Glomeromycotina 做了有意加权。")
    _pr = [r for r in rows if r["group"] == "优先"]
    _gcf = sum(1 for r in _pr if r["acc"].startswith("GCF_"))
    _jgi = sum(1 for r in _pr if ("Joint Genome Institute" in r["src"]) or r["acc"].startswith("JGI:"))
    _rna = sum(1 for r in _pr if not r["rna"].startswith("无"))
    _src = {}
    for r in _pr:
        k = ("JGI" if (("Joint Genome Institute" in r["src"]) or r["acc"].startswith("JGI:")) else
             "Kew (RBG Kew)" if "Kew" in r["src"] else
             "Broad Institute" if "Broad" in r["src"] else
             "NHGRI" if "NHGRI" in r["src"] else "NIBB" if "NIBB" in r["src"] else r["src"])
        _src[k] = _src.get(k, 0) + 1
    bullet(f"优先组 10 株的来源构成：RefSeq（GCF_）参考基因组 {_gcf} 株；"
           f"既有注释来自 JGI 的 {_jgi} 株（含 1 株 MycoCosm 门户直取），其余为 "
           + "、".join(f"{k} {v} 株" for k, v in _src.items() if k != "JGI") + "。")
    _norun = "、".join(r["name"] for r in _pr if r["rna"].startswith("无"))
    bullet(f"优先组 {_rna} 株在队列内保有 RNA-seq run 可直接调用（BRAKER3 侧可直接复用）；"
           f"无保留 run 的 {10 - _rna} 株（{_norun}）需在 BRAKER3 与 BRAKER2 两层之间对齐后再比较。")
    bullet("「既有注释有 RNA-seq 支持」这一条比预想的稀缺，且分层很细。已确证有 RNA-seq 证据的："
           "Mucor velutinosus NIH1002（16 个生长条件）、Rhizophagus irregularis（短读 + 长读双证据）、"
           "Rhizophagus clarus（BRAKER v2.1.2 + RNA-seq，但所用数据集未点名）。"
           "JGI 各株的流程以转录组数据为标准输入（BLAT + GMAP 映射后驱动 FGENESH / EST_map / combest / estExt），"
           "但逐株是否采用未在公开元数据中标注。"
           "已确证没有的：Lichtheimia ornata（BRAKER2 蛋白模式）、EBP-Nor 各株（纯蛋白同源）、"
           "Mucor flavus / Helicostylum pulchrum（Funannotate 1.8.15 默认）、Linnemannia hyalina（Funannotate 1.8.1 全默认）、"
           "Gigaspora 两株（GeneMark-ES 纯从头预测）。完整分层见第 1.5 节，这直接决定了后续必须按证据类型分层报告。")
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
           ["⑥", "按亚门配额与来源可信度人工定序 → 优先 10 + 备选 10", 20]],
          widths=[1.6, 11.4, 3.6], font=9)
    para(f"放宽档（BUSCO C ≥ 98.0% 且 N50 ≥ 500 kb 或染色体级）可得到 {stats['n_loose']} 株，"
         f"即本报告的扩展替补池；若优先组在后续核查中需要替换，从这里按亚门就近取用。", size=9)

    doc.add_heading("1.5　既有注释的证据类型分层（本报告最重要的前置结论）", level=3)
    para("这一节回答一个比「哪株基因组好」更要命的问题：作为基准的那份注释，本身是用什么证据、由什么流程产生的。"
         "结论直接影响后续所有指标的解释方式。以下每一条都核实到一手来源（论文原文、ENA/NCBI 记录、流程仓库），"
         "无法核实的项如实标注「未验证」。")
    table(["物种", "Accession", "注释流程", "RNA-seq 证据", "人工审校", "作为基准的可信度", "来源"],
          [list(r) for r in EVIDENCE_TIERS],
          widths=[3.0, 2.5, 3.0, 4.3, 2.2, 2.2, 3.0], font=7.5)
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
           "并用 BUSCO、OMArk、Swiss-Prot/InterPro 支持率等不依赖参考基因集的指标作为主证据。详见第 6.10 节。")

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
    para("优先组是全套比对的主战场：建议把批次一任务 1.3–1.5 的全部指标先在这 10 株上跑完，"
         "确认流程与判据稳定后，再决定是否扩到备选组。")
    table(["编号", "物种", "菌株", "亚门", "Accession", "装配", "N50(kb)", "BUSCO C(%)",
           "既有注释来源", "基因数"],
          [[r["no"], r["name"], r["strain"], r["sub_short"], r["acc"], r["level"], r["n50"],
            r["busco"], r["src_short"], r["genes"]] for r in rows if r["group"] == "优先"],
          widths=[1.1, 3.5, 1.9, 2.0, 2.7, 1.5, 1.2, 1.3, 3.0, 1.3], font=7.5)
    for r in rows:
        if r["group"] == "优先":
            species_block(r["tag"], r["name"], r["strain"], r["sub"], r["acc"], r["src"], r["level"],
                          r["size"], r["n50"], r["busco"], r["genes"], r["rna"], r["orthodb"],
                          r["isoform"], r["reason"], r["risk"])

    doc.add_heading("4　备选组（10 株）", level=2)
    para("备选组的用途有三个：一是优先组某一株在任务 1.2/1.3 暴露问题时的就近替换；"
         "二是把'基因组好但注释来源是自动流程'的样本单列，用于分离装配质量与注释质量两个变量；"
         "三是把亚门内多样性再扩一层，避免所有结论都建立在少数几个属上。")
    table(["编号", "物种", "菌株", "亚门", "Accession", "装配", "N50(kb)", "BUSCO C(%)",
           "既有注释来源", "基因数"],
          [[r["no"], r["name"], r["strain"], r["sub_short"], r["acc"], r["level"], r["n50"],
            r["busco"], r["src_short"], r["genes"]] for r in rows if r["group"] == "备选"],
          widths=[1.1, 3.5, 1.9, 2.0, 2.7, 1.5, 1.2, 1.3, 3.0, 1.3], font=7.5)
    for r in rows:
        if r["group"] == "备选":
            species_block(r["tag"], r["name"], r["strain"], r["sub"], r["acc"], r["src"], r["level"],
                          r["size"], r["n50"], r["busco"], r["genes"], r["rna"], r["orthodb"],
                          r["isoform"], r["reason"], r["risk"])

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
            "20 株的 accession、来源、基因数、基因密度",
            "五项指标：精确一致率、双向 RBH、互命中对一致度、拆分/合并、各自独有蛋白"],
           ["1.6 结果整理与呈现",
            "本报告的亚门与来源分组口径",
            "Table 5 按亚门与证据库归属分组；措辞统一用'一致性/吻合程度'，不用'准确性'"]],
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

def main():
    m = load_master()
    membership = orthodb_members()
    ncbi = fetch_ncbi([r[0] for r in PRIORITY + BACKUP] + [e[0] for e in EXTENDED])
    m["busco"] = pd.to_numeric(m["BUSCO C_muco"], errors="coerce")
    m["n50kb"] = pd.to_numeric(m["assembly_N50_kb"], errors="coerce").astype(float)
    # 参考注释异构体负荷（由已取回的参考 GFF3 统计）
    iso_path = os.path.join(OUTDIR, "reference_gff_isoform_stats.json")
    pstat_path = os.path.join(OUTDIR, "reference_protein_stats.json")
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
            risk=risk or "无特别风险；按亚门分组报告即可。",
        )

    rows = []
    for i, (acc, tag, reason, risk) in enumerate(PRIORITY, 1):
        rows.append(pack(acc, "优先", f"P{i}", reason, risk))
    for i, (acc, sub, reason) in enumerate(BACKUP, 1):
        rows.append(pack(acc, "备选", f"B{i}", reason, ""))

    stats = {
        "n_total": len(m),
        "n_has_ann": int((m["原始蛋白注释"] == "有").sum()),
        "n_attrib": int(((m["原始蛋白注释"] == "有") &
                         (m["来源库"].eq("JGI") | m["主标识"].astype(str).str.startswith(("GCA_", "GCF_")))).sum()),
        "n_busco985": int(((m["原始蛋白注释"] == "有") & (m["busco"] >= 98.5)).sum()),
        "n_final": int(((m["原始蛋白注释"] == "有") & (m["busco"] >= 98.5) &
                        ((m["n50kb"] >= 800) | (m["assembly 水平（来源词）"].isin(["Chromosome", "Complete Genome"])))).sum()),
        "n_loose": int(((m["原始蛋白注释"] == "有") & (m["busco"] >= 98.0) &
                        ((m["n50kb"] >= 500) | (m["assembly 水平（来源词）"].isin(["Chromosome", "Complete Genome"])))).sum()),
    }

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
    csv_path = os.path.join(OUTDIR, "基准物种20株_名单.csv")
    with open(csv_path, "w", encoding="utf-8-sig") as fh:
        fh.write("组别,标签,No,物种,菌株,亚门,Assembly_accession,装配水平,基因组大小_Mb,N50_kb,"
                 "基因组BUSCO_C_pct,既有注释提交方,既有注释基因数,OrthoDB_v12.2_证据库,参考注释异构体,队列内RNA-seq\n")
        for r in rows:
            fh.write(",".join('"%s"' % str(r[k]).replace('"', "'") for k in
                              ["group", "tag", "no", "name", "strain", "sub", "acc", "level", "size",
                               "n50", "busco", "src", "genes", "orthodb", "isoform", "rna"]) + "\n")
    print("written:", csv_path)

    notes = json.load(open(os.path.join(OUTDIR, "risk_notes.json")))
    path = build_docx(rows, stats, [(n["heading"], n["bullets"]) for n in notes], REF_CHECK)
    print("written:", path)
    print("stats:", stats)


if __name__ == "__main__":
    main()
