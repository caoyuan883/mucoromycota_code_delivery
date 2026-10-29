#!/bin/bash
# ============================================================================
# Step 13 — 把可丢弃的中间产物搬进 _TRASH_pending（不删除！）
#   DRY_RUN=1 只统计将搬走多少
#   同一文件系统内 mv 是瞬时的，随时可以搬回来
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

DRY_RUN=${DRY_RUN:-0}
TRASH=$BATCH_ROOT/_TRASH_pending
LOG=$BATCH_ROOT/_cleanup_stage.log
mkdir -p $TRASH

# 整目录搬走
DROP_DIRS="Spaln GeneMark-ES genome_split bbc.backup_pre_merge mb_downloads"
# GeneMark-EP 内部子目录搬走（保留最终 gtf）
DROP_GM_SUB="run data output info"
# 单文件搬走
DROP_FILES="nuc.fasta prevHints.gff augustus.hints_iter1.gff augustus.hints.gff augustus.hints.tmp.gtf braker.gtf_temp braker.aa.backup_pre_merge braker.gtf.backup_pre_merge genome.fa"

if [[ "$DRY_RUN" == "1" ]]; then
    echo "===== DRY RUN — 只统计将搬走的体积 ====="
    tot=0
    for d in $BATCH_ROOT/*/; do
        sp=$(basename $d); [[ "$sp" == _* ]] && continue
        for x in $DROP_DIRS;  do [ -e "$d/$x" ] && tot=$((tot + $(du -sb "$d/$x" 2>/dev/null | cut -f1))); done
        for x in $DROP_GM_SUB; do [ -e "$d/GeneMark-EP/$x" ] && tot=$((tot + $(du -sb "$d/GeneMark-EP/$x" 2>/dev/null | cut -f1))); done
        for x in $DROP_FILES; do [ -e "$d/$x" ] && tot=$((tot + $(stat -c%s "$d/$x" 2>/dev/null || echo 0))); done
    done
    echo "  将搬入待删区: $(numfmt --to=iec $tot)"
    echo "  当前批次目录: $(du -sh --exclude=_TRASH_pending $BATCH_ROOT 2>/dev/null | cut -f1)"
    echo "  当前待删区:   $(du -sh $TRASH 2>/dev/null | cut -f1)"
    exit 0
fi

echo "[$(ts)] ===== 搬移可丢弃中间产物 =====" | tee -a $LOG
n=0
for d in $BATCH_ROOT/*/; do
    sp=$(basename $d); [[ "$sp" == _* ]] && continue
    n=$((n+1))
    T=$TRASH/$sp
    mkdir -p $T
    for x in $DROP_DIRS;  do [ -e "$d/$x" ] && mv "$d/$x" "$T/" 2>/dev/null; done
    if [ -d "$d/GeneMark-EP" ]; then
        mkdir -p $T/GeneMark-EP
        for x in $DROP_GM_SUB; do [ -e "$d/GeneMark-EP/$x" ] && mv "$d/GeneMark-EP/$x" "$T/GeneMark-EP/" 2>/dev/null; done
        rmdir $T/GeneMark-EP 2>/dev/null
    fi
    for x in $DROP_FILES; do [ -e "$d/$x" ] && mv "$d/$x" "$T/" 2>/dev/null; done
    rmdir $T 2>/dev/null
    echo "[$(ts)] [$n/127] $sp" >> $LOG
done

# 批次级冗余（已并入 _archive_intermediate/_batch_level.tar.gz）
for x in _bbc_all _bbc_merge _packages; do
    [ -e "$BATCH_ROOT/$x" ] && mv "$BATCH_ROOT/$x" "$TRASH/" 2>/dev/null && echo "[$(ts)] 批次级: $x" >> $LOG
done

echo "" | tee -a $LOG
echo "[$(ts)] ===== 搬移完成 =====" | tee -a $LOG
echo "  处理物种: $n" | tee -a $LOG
echo "  待删区:   $(du -sh $TRASH 2>/dev/null | cut -f1)" | tee -a $LOG
echo "  批次目录: $(du -sh --exclude=_TRASH_pending $BATCH_ROOT 2>/dev/null | cut -f1)" | tee -a $LOG
echo "  可用空间: $(df -h /media/desk16 | awk 'NR==2{print $4}')" | tee -a $LOG
