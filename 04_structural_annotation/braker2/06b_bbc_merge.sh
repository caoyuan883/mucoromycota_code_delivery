#!/bin/bash
# ============================================================================
# Step 06b v3 — Missing >5% 物种跑 BRAKER 内部 best_by_compleasm.py
# v3 修复:
#   1. 固定 cwd 到 mb_downloads 的父目录
#   2. 每株跑前保障 placement_files.done / lineage .done 标记（绕过联网）
#   3. 清 placement_files.tmp 僵尸锁
#   4. 去掉 -b（回到与冒烟成功一致的参数）
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

# [release refactor] historical: tool paths under /media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src
BBC_SCRIPT=$ENV_ROOT/braker2_bundle_transfer/src/TSEBRA/bin/best_by_compleasm.py
GETANNO=$ENV_ROOT/braker2_bundle_transfer/src/AUGUSTUS_src/scripts/getAnnoFastaFromJoingenes.py
TSEBRA=$ENV_ROOT/braker2_bundle_transfer/src/TSEBRA/bin/tsebra.py

# 关键：compleasm 用相对路径 mb_downloads，所以固定 cwd
# [release refactor] historical: CWD_FIX=/media/desk16/tly8529/braker2_batch_v3 (deployed copy of this folder; compleasm uses ./mb_downloads)
CWD_FIX=${BBC_CWD:-$here}
MB=$CWD_FIX/mb_downloads

MERGE_ROOT=$BATCH_ROOT/_bbc_merge
MERGE_LOG=$BATCH_ROOT/_bbc_merge.log
MERGE_MANIFEST=$BATCH_ROOT/_bbc_merge_manifest.tsv

mkdir -p $MERGE_ROOT $MB
: > $MERGE_LOG

log_m() { echo "[$(ts)] [MERGE] $*" | tee -a $MERGE_LOG; }

# ---- 保障 compleasm 不联网 -----------------------------------------------
ensure_offline() {
    \rm -f $MB/placement_files.tmp
    mkdir -p $MB/placement_files
    touch $MB/placement_files.done
    touch $MB/eukaryota_odb12.done
    touch $MB/mucoromycota_odb12.done
    touch $MB/file_versions.tsv.done
}

MERGE_LIST=$MERGE_ROOT/_species_to_merge.txt
awk -F'\t' 'NR>1 && $9!="NA" && $9+0 > 5 {print $1}' $BATCH_ROOT/_manifest.tsv > $MERGE_LIST
n_todo=$(wc -l < $MERGE_LIST)

log_m "===== BRAKER2 v3 best_by_compleasm 合并 v3 ====="
log_m "  待处理物种数: $n_todo"
log_m "  固定 cwd:     $CWD_FIX"
log_m "  mb_downloads: $MB"
log_m "  预计每株:     20-30 min"
log_m "  预计总耗时:   $((n_todo * 25)) min"
log_m ""
log_m "===== 物种列表 ====="
cat $MERGE_LIST | tee -a $MERGE_LOG

echo -e "species\tstatus\tbraker_miss\tgenemark_miss\taugustus_miss\tfinal_gtf\tfinal_prot\tdecision\telapsed_s" > $MERGE_MANIFEST

cd $CWD_FIX
ensure_offline

merge_one() {
    local sp=$1
    local idx=$2
    local total=$3

    local wd=$BATCH_ROOT/$sp
    local genome=$wd/genome.fa
    local tmp=$MERGE_ROOT/$sp
    local log=$MERGE_ROOT/${sp}.log

    if [[ ! -s "$genome" ]]; then
        log_m "  [$idx/$total] [$sp] ❌ 缺 genome.fa"
        echo -e "$sp\tSKIP_NO_GENOME\tNA\tNA\tNA\tno\tNA\tNA\t0" >> $MERGE_MANIFEST
        return 1
    fi
    if [[ ! -s "$wd/braker.aa" ]]; then
        log_m "  [$idx/$total] [$sp] ❌ 缺 braker.aa"
        echo -e "$sp\tSKIP_NO_BRAKER\tNA\tNA\tNA\tno\tNA\tNA\t0" >> $MERGE_MANIFEST
        return 1
    fi
    if [[ ! -s "$wd/GeneMark-EP/genemark.gtf" ]]; then
        log_m "  [$idx/$total] [$sp] ❌ 缺 GeneMark-EP/genemark.gtf"
        echo -e "$sp\tSKIP_NO_GENEMARK\tNA\tNA\tNA\tno\tNA\tNA\t0" >> $MERGE_MANIFEST
        return 1
    fi

    # 备份（只备份一次，绝不删除原件）
    [[ ! -d "$wd/bbc.backup_pre_merge" ]] && cp -r $wd/bbc $wd/bbc.backup_pre_merge 2>/dev/null
    [[ ! -f "$wd/braker.aa.backup_pre_merge" ]] && cp $wd/braker.aa $wd/braker.aa.backup_pre_merge
    [[ ! -f "$wd/braker.gtf.backup_pre_merge" ]] && cp $wd/braker.gtf $wd/braker.gtf.backup_pre_merge

    \rm -rf $tmp
    mkdir -p $tmp

    ensure_offline

    log_m "  [$idx/$total] [$sp] 启动合并 (预计 25 min)"

    local t0=$(date +%s)
    python3 $BBC_SCRIPT \
        -m $tmp \
        -c $COMPLEASM_PATH/compleasm.py \
        -d $wd \
        -y $TSEBRA \
        -f $GETANNO \
        -g $genome \
        -t 8 \
        -p $BUSCO_LINEAGE_DIR \
        > $log 2>&1
    local rc=$?
    local t1=$(date +%s)
    local elapsed=$((t1 - t0))

    ensure_offline

    if [[ "$rc" -ne 0 ]]; then
        log_m "  [$idx/$total] [$sp] ❌ 失败 rc=$rc elapsed=${elapsed}s"
        log_m "    错误尾部:"
        tail -5 $log | sed 's/^/      /' | tee -a $MERGE_LOG
        echo -e "$sp\tFAIL_rc$rc\tNA\tNA\tNA\tno\tNA\tNA\t$elapsed" >> $MERGE_MANIFEST
        return 1
    fi

    local br_miss=$(grep "BRAKER is missing" $log | head -1 | awk '{print $4}')
    local gm_miss=$(grep "GeneMark is missing" $log | head -1 | awk '{print $4}')
    local au_miss=$(grep "Augustus is missing" $log | head -1 | awk '{print $4}')
    [[ -z "$br_miss" ]] && br_miss="NA"
    [[ -z "$gm_miss" ]] && gm_miss="NA"
    [[ -z "$au_miss" ]] && au_miss="NA"

    local decision="unknown"
    if grep -q "is the best one" $log; then
        decision="keep_braker_best"
    elif grep -q "there are no BUSCOs to be added" $log; then
        decision="keep_braker_no_add"
    elif grep -q "will be kept" $log; then
        decision="keep_braker"
    elif grep -qi "merge\|added to the braker" $log; then
        decision="merged"
    fi

    local final_gen="no"
    local final_n="NA"
    if [[ -s "$wd/braker_final.gtf" ]]; then
        final_gen="yes"
        [[ -s "$wd/braker_final.aa" ]] && final_n=$(grep -c '^>' $wd/braker_final.aa)
    fi

    log_m "  [$idx/$total] [$sp] ✅ br=${br_miss}% gm=${gm_miss}% au=${au_miss}% decision=$decision final=$final_gen n=$final_n ${elapsed}s"
    echo -e "$sp\tDONE\t$br_miss\t$gm_miss\t$au_miss\t$final_gen\t$final_n\t$decision\t$elapsed" >> $MERGE_MANIFEST
    return 0
}

i=0
for sp in $(cat $MERGE_LIST); do
    i=$((i + 1))
    merge_one "$sp" "$i" "$n_todo"
done

log_m ""
log_m "===== 合并完成 ====="
n_done=$(awk -F'\t' 'NR>1 && $2=="DONE"' $MERGE_MANIFEST | wc -l)
n_fail=$(awk -F'\t' 'NR>1 && $2 ~ /^FAIL/' $MERGE_MANIFEST | wc -l)
n_skip=$(awk -F'\t' 'NR>1 && $2 ~ /^SKIP/' $MERGE_MANIFEST | wc -l)
n_final=$(awk -F'\t' 'NR>1 && $6=="yes"' $MERGE_MANIFEST | wc -l)
log_m "  DONE:                  $n_done"
log_m "  FAIL:                  $n_fail"
log_m "  SKIP:                  $n_skip"
log_m "  生成 braker_final.gtf: $n_final"
log_m ""
column -t -s $'\t' $MERGE_MANIFEST | tee -a $MERGE_LOG
