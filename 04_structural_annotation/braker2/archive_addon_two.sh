#!/bin/bash
# ============================================================================
# Step 16 — 只归档新增 2 株的中间产物（保留清单与 archive_intermediates.sh 完全一致）
#   输出: $BATCH_ROOT/_archive_intermediate_2_new/
#   额外: _batch_level_update.tar.gz —— 更新后的 manifest/日志/脚本 + 新 2 株的 _bbc_full
#   安全: 只读源目录，只写新目录，不碰 _archive_intermediate
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

SPECIES="Cryptocalciella_humilis Rhizopus_sp._TYR1"
ARCH=$BATCH_ROOT/_archive_intermediate_2_new
LOG=$BATCH_ROOT/_archive_2_new.log
IDX=$ARCH/index.tsv

mkdir -p $ARCH
log_a() { echo "[$(ts)] [ARCH2] $*" | tee -a $LOG; }

# 与 Step 12 完全相同的保留清单
KEEP_DIRS="species errors"
KEEP_FILES="prothint.gff evidence.gff prothint_augustus.gff augustus.hints.gtf augustus.hints.aa genome_header.map what-to-cite.txt"
KEEP_GM="genemark.gtf genemark.f.good.gtf"
KEEP_BBC="genemark.aa"

echo -e "species\tarchive\tsize\tmd5" > $IDX

log_a "===== 新增 2 株中间产物归档 ====="
log_a "  输出: $ARCH"

ok=0
for sp in $SPECIES; do
    d=$BATCH_ROOT/$sp
    [ -d "$d" ] || { log_a "  [$sp] 工作目录不存在，跳过"; continue; }
    stage=$ARCH/$sp
    \rm -rf $stage
    mkdir -p $stage

    for x in $KEEP_DIRS;  do [ -e "$d/$x" ] && cp -r "$d/$x" "$stage/"; done
    for x in $KEEP_FILES; do [ -e "$d/$x" ] && cp "$d/$x" "$stage/"; done
    if [ -d "$d/GeneMark-EP" ]; then
        mkdir -p $stage/GeneMark-EP
        for x in $KEEP_GM gmes.log run.cfg; do
            [ -e "$d/GeneMark-EP/$x" ] && cp "$d/GeneMark-EP/$x" "$stage/GeneMark-EP/"
        done
    fi
    if [ -d "$d/bbc" ]; then
        mkdir -p $stage/bbc
        for x in $KEEP_BBC; do [ -e "$d/bbc/$x" ] && cp "$d/bbc/$x" "$stage/bbc/"; done
        for s in braker augustus genemark; do
            [ -d "$d/bbc/$s" ] && { mkdir -p $stage/bbc/$s; cp $d/bbc/$s/*.txt $d/bbc/$s/*.tsv $stage/bbc/$s/ 2>/dev/null; }
        done
    fi
    for f in $d/*.stdout $d/*.stderr; do [ -f "$f" ] && cp "$f" "$stage/"; done

    # 后置择优的完整结果（06c）
    [ -d "$BATCH_ROOT/_bbc_full/$sp" ] && cp -r "$BATCH_ROOT/_bbc_full/$sp" "$stage/bbc_full_06c"
    [ -f "$BATCH_ROOT/_bbc_full/${sp}.log" ] && cp "$BATCH_ROOT/_bbc_full/${sp}.log" "$stage/bbc_full_06c.log"

    cat > $stage/ARCHIVE_NOTE.txt <<EOF
Intermediate files retained for ${sp//_/ }.

This genome was added to the dataset in August 2026, after the initial batch of
125 genomes. It was processed with an identical BRAKER2 command line, the same
protein evidence set (Fungi.clean.fa) and the same BUSCO lineage
(mucoromycota_odb12), so the annotation is directly comparable to the rest of
the resource.

Retained
  species/              AUGUSTUS species parameters trained for this genome
  GeneMark-EP/          final GeneMark-EP+ gene models and run log
  augustus.hints.*      AUGUSTUS candidate gene set used in gene-set selection
  prothint.gff          ProtHint protein evidence
  bbc/                  compleasm results written by the BRAKER run
  bbc_full_06c/         compleasm results of the post-hoc gene-set selection
  errors/, *.stdout/err complete runtime records

Not retained (regenerable or superseded)
  proteins.fa, diamond/ identical copy of the ProtHint input database
                        (Fungi.clean.fa) and its DIAMOND index
  Spaln/                Spaln alignment scratch
  GeneMark-ES/          superseded by GeneMark-EP+
  GeneMark-EP/{run,data,output}   GeneMark scratch directories
  genome_split/, nuc.fasta, prevHints.gff, *_iter1.gff, *.tmp.gtf
                        BRAKER iteration scratch
  genome.fa             soft-masked assembly, kept in 03_softmasked_genomes_v3
  *.backup_pre_merge    duplicates of files in the core release package

Final annotation products are distributed in ${sp}_braker2_v4.tar.gz.
EOF

    tar -czf $ARCH/${sp}_intermediate.tar.gz -C $ARCH $sp 2>/dev/null && \rm -rf $stage
    if [ -s "$ARCH/${sp}_intermediate.tar.gz" ]; then
        sz=$(du -h $ARCH/${sp}_intermediate.tar.gz | cut -f1)
        m5=$(md5sum $ARCH/${sp}_intermediate.tar.gz | cut -d' ' -f1)
        echo -e "$sp\t${sp}_intermediate.tar.gz\t$sz\t$m5" >> $IDX
        ok=$((ok+1))
        log_a "  [$sp] OK $sz"
    else
        log_a "  [$sp] FAIL"
        \rm -rf $stage
    fi
done

# ---- 更新后的批次级材料（manifest 已从 125 变 127）----
log_a "打包更新后的批次级材料 ..."
BST=$ARCH/_batch_level_update
\rm -rf $BST; mkdir -p $BST
cp $BATCH_ROOT/_manifest_final.tsv $BATCH_ROOT/_manifest_final_report.txt $BST/ 2>/dev/null
cp $BATCH_ROOT/_manifest_v4.tsv $BATCH_ROOT/_manifest_v4_report.txt $BST/ 2>/dev/null
cp $BATCH_ROOT/_busco_consistency.tsv $BATCH_ROOT/_bbc_full_manifest.tsv $BST/ 2>/dev/null
cp $BATCH_ROOT/_packages_v4_index.tsv $BATCH_ROOT/_archive_index.tsv $BST/ 2>/dev/null
cp $BATCH_ROOT/_bbc_full.log $BATCH_ROOT/_package_v4.log $BATCH_ROOT/_gff3_backfill.log $BST/ 2>/dev/null
cp $BATCH_ROOT/_addon_260807.log $BATCH_ROOT/_archive_2_new.log $BST/ 2>/dev/null
cp $BATCH_ROOT/_logs/Cryptocalciella_humilis.* $BATCH_ROOT/_logs/Rhizopus_sp._TYR1.* $BST/ 2>/dev/null
cp -r $here $BST/pipeline_scripts 2>/dev/null
\rm -f $BST/pipeline_scripts/*.bak_125 2>/dev/null
tar -czf $ARCH/_batch_level_update.tar.gz -C $ARCH _batch_level_update 2>/dev/null && \rm -rf $BST

( cd $ARCH && md5sum *.tar.gz > MD5SUMS.txt )

log_a ""
log_a "===== 完成 ====="
log_a "  物种包: $ok / 2"
log_a "  总大小: $(du -sh $ARCH | cut -f1)"
ls -lh $ARCH/*.tar.gz | awk '{print "    "$5"  "$9}' | tee -a $LOG
