#!/bin/bash
# ============================================================================
# Step 12 — 精选中间产物归档
#  设 DRY_RUN=1 只统计体积，不生成任何文件
#  保留原则: 可复现所必需 + 难以再生 + 全部运行记录
#  安全: 只读源目录，只写 _archive_intermediate/
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

DRY_RUN=${DRY_RUN:-0}
ARCH=$BATCH_ROOT/_archive_intermediate
LOG=$BATCH_ROOT/_archive.log
IDX=$BATCH_ROOT/_archive_index.tsv

mkdir -p $ARCH
log_a() { echo "[$(ts)] [ARCH] $*" | tee -a $LOG; }

# 每株保留清单
KEEP_DIRS="species errors"
KEEP_FILES="prothint.gff evidence.gff prothint_augustus.gff augustus.hints.gtf augustus.hints.aa genome_header.map what-to-cite.txt"
KEEP_GM="genemark.gtf genemark.f.good.gtf"
KEEP_BBC="genemark.aa"

if [[ "$DRY_RUN" == "1" ]]; then
    echo "===== DRY RUN — 只统计，不写文件 ====="
    total=0
    for d in $BATCH_ROOT/*/; do
        sp=$(basename $d); [[ "$sp" == _* ]] && continue
        s=0
        for x in $KEEP_DIRS; do [ -e "$d/$x" ] && s=$((s + $(du -sb "$d/$x" 2>/dev/null | cut -f1))); done
        for x in $KEEP_FILES; do [ -e "$d/$x" ] && s=$((s + $(stat -c%s "$d/$x" 2>/dev/null || echo 0))); done
        for x in $KEEP_GM;  do [ -e "$d/GeneMark-EP/$x" ] && s=$((s + $(stat -c%s "$d/GeneMark-EP/$x"))); done
        for x in $KEEP_BBC; do [ -e "$d/bbc/$x" ] && s=$((s + $(stat -c%s "$d/bbc/$x"))); done
        for f in $d/*.stdout $d/*.stderr; do [ -f "$f" ] && s=$((s + $(stat -c%s "$f"))); done
        total=$((total + s))
    done
    echo "  每株保留合计（未压缩）: $(numfmt --to=iec $total)"
    echo "  _bbc_full:              $(du -sh $BATCH_ROOT/_bbc_full 2>/dev/null | cut -f1)"
    echo "  批次日志清单:            $(du -ch $BATCH_ROOT/_*.tsv $BATCH_ROOT/_*.log $BATCH_ROOT/_*.txt 2>/dev/null | tail -1 | cut -f1)"
    echo ""
    echo "  预计压缩后约为未压缩的 30-40%"
    exit 0
fi

echo -e "species\tarchive\tsize\tmd5" > $IDX

log_a "===== 精选中间产物归档 ====="
log_a "  输出: $ARCH"

n=0; ok=0
for d in $BATCH_ROOT/*/; do
    sp=$(basename $d); [[ "$sp" == _* ]] && continue
    n=$((n+1))
    stage=$ARCH/$sp
    \rm -rf $stage
    mkdir -p $stage

    for x in $KEEP_DIRS;  do [ -e "$d/$x" ] && cp -r "$d/$x" "$stage/"; done
    for x in $KEEP_FILES; do [ -e "$d/$x" ] && cp "$d/$x" "$stage/"; done
    if [ -d "$d/GeneMark-EP" ]; then
        mkdir -p $stage/GeneMark-EP
        for x in $KEEP_GM; do [ -e "$d/GeneMark-EP/$x" ] && cp "$d/GeneMark-EP/$x" "$stage/GeneMark-EP/"; done
        for x in gmes.log run.cfg; do [ -e "$d/GeneMark-EP/$x" ] && cp "$d/GeneMark-EP/$x" "$stage/GeneMark-EP/"; done
    fi
    if [ -d "$d/bbc" ]; then
        mkdir -p $stage/bbc
        for x in $KEEP_BBC; do [ -e "$d/bbc/$x" ] && cp "$d/bbc/$x" "$stage/bbc/"; done
        for s in braker augustus genemark; do
            [ -d "$d/bbc/$s" ] && { mkdir -p $stage/bbc/$s; cp $d/bbc/$s/*.txt $d/bbc/$s/*.tsv $stage/bbc/$s/ 2>/dev/null; }
        done
    fi
    for f in $d/*.stdout $d/*.stderr; do [ -f "$f" ] && cp "$f" "$stage/"; done

    cat > $stage/ARCHIVE_NOTE.txt <<EOF
Intermediate files retained for ${sp//_/ }.

Retained
  species/              AUGUSTUS species parameters trained for this genome
  GeneMark-EP/          final GeneMark-EP+ gene models and run log
  augustus.hints.*      AUGUSTUS candidate gene set used in gene-set selection
  prothint.gff          ProtHint protein evidence
  bbc/                  compleasm results of the gene-set selection step
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
        log_a "  [$n/127] [$sp] OK $sz"
    else
        log_a "  [$n/127] [$sp] FAIL"
        \rm -rf $stage
    fi
done

# ---- 批次级归档 ----
log_a "打包批次级材料 ..."
BSTAGE=$ARCH/_batch_level
\rm -rf $BSTAGE; mkdir -p $BSTAGE
cp -r $BATCH_ROOT/_bbc_full $BSTAGE/ 2>/dev/null
cp $BATCH_ROOT/_*.tsv $BATCH_ROOT/_*.log $BATCH_ROOT/_*.txt $BSTAGE/ 2>/dev/null
cp -r $here $BSTAGE/pipeline_scripts 2>/dev/null
tar -czf $ARCH/_batch_level.tar.gz -C $ARCH _batch_level 2>/dev/null && \rm -rf $BSTAGE

( cd $ARCH && md5sum *.tar.gz > MD5SUMS.txt )

log_a ""
log_a "===== 归档完成 ====="
log_a "  物种包: $ok / $n"
log_a "  总大小: $(du -sh $ARCH | cut -f1)"
log_a "  索引:   $IDX"
