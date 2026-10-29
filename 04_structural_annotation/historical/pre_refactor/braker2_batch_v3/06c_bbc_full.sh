#!/bin/bash
# ============================================================================
# Step 06c — 对全部 125 株跑 BRAKER 内部 best_by_compleasm.py
# 安全原则: 只写新目录 _bbc_full；断点续跑；串行
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

BBC_SCRIPT=/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/TSEBRA/bin/best_by_compleasm.py
GETANNO=/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/AUGUSTUS_src/scripts/getAnnoFastaFromJoingenes.py
TSEBRA=/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/TSEBRA/bin/tsebra.py
BUSCO_LIB=/media/desk16/tly8529/mucoromycota_odb12

CWD_FIX=/media/desk16/tly8529/braker2_batch_v3
MB=$CWD_FIX/mb_downloads

FULL_ROOT=$BATCH_ROOT/_bbc_full
STATE_DIR=$FULL_ROOT/_state
FULL_LOG=$BATCH_ROOT/_bbc_full.log
FULL_MANIFEST=$BATCH_ROOT/_bbc_full_manifest.tsv

mkdir -p $FULL_ROOT $STATE_DIR $MB

log_f() { echo "[$(ts)] [FULL] $*" | tee -a $FULL_LOG; }

ensure_offline() {
    \rm -f $MB/placement_files.tmp
    mkdir -p $MB/placement_files
    touch $MB/placement_files.done
    touch $MB/eukaryota_odb12.done
    touch $MB/mucoromycota_odb12.done
    touch $MB/file_versions.tsv.done
}

getm() { grep '^M:' "$1" 2>/dev/null | head -1 | sed 's/^M://; s/%.*//; s/^ *//'; }
gets() { grep '^S:' "$1" 2>/dev/null | head -1 | sed 's/^S://; s/%.*//; s/^ *//'; }
getd() { grep '^D:' "$1" 2>/dev/null | head -1 | sed 's/^D://; s/%.*//; s/^ *//'; }
getf() { grep '^F:' "$1" 2>/dev/null | head -1 | sed 's/^F://; s/%.*//; s/^ *//'; }

n_total=$(wc -l < $SPECIES_LIST_FILE)

if [[ ! -s "$FULL_MANIFEST" ]]; then
  echo -e "species\tstatus\tmode\tbraker_M\taugustus_M\tgenemark_M\tbraker_S\tbraker_D\tbraker_F\tfinal_gtf\tfinal_prot\tdecision\telapsed_s" > $FULL_MANIFEST
fi

log_f "===== BRAKER2 v3 全量 best_by_compleasm (06c) ====="
log_f "  物种总数:   $n_total"
log_f "  输出目录:   $FULL_ROOT"
log_f "  BUSCO 库:   $BUSCO_LIB"
log_f "  固定 cwd:   $CWD_FIX"
log_f "  预计耗时:   约 $((n_total * 11 / 60)) 小时"

cd $CWD_FIX
ensure_offline

run_one() {
    local sp=$1 idx=$2
    local wd=$BATCH_ROOT/$sp
    local genome=$wd/genome.fa
    local tmp=$FULL_ROOT/$sp
    local log=$FULL_ROOT/${sp}.log
    local st=$STATE_DIR/${sp}.state

    if [[ -f "$st" ]] && grep -q '^DONE' "$st"; then
        log_f "  [$idx/$n_total] [$sp] 已完成，跳过"
        return 0
    fi

    if [[ ! -s "$genome" ]]; then
        log_f "  [$idx/$n_total] [$sp] SKIP 缺 genome.fa"
        echo "SKIP_NO_GENOME" > $st
        echo -e "$sp\tSKIP_NO_GENOME\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno\tNA\tNA\t0" >> $FULL_MANIFEST
        return 1
    fi
    if [[ ! -s "$wd/braker.aa" ]]; then
        log_f "  [$idx/$n_total] [$sp] SKIP 缺 braker.aa"
        echo "SKIP_NO_BRAKER" > $st
        echo -e "$sp\tSKIP_NO_BRAKER\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno\tNA\tNA\t0" >> $FULL_MANIFEST
        return 1
    fi

    [[ -d "$wd/bbc" && ! -d "$wd/bbc.backup_pre_merge" ]] && cp -r $wd/bbc $wd/bbc.backup_pre_merge 2>/dev/null
    [[ ! -f "$wd/braker.aa.backup_pre_merge" ]] && cp $wd/braker.aa $wd/braker.aa.backup_pre_merge
    [[ ! -f "$wd/braker.gtf.backup_pre_merge" ]] && cp $wd/braker.gtf $wd/braker.gtf.backup_pre_merge

    \rm -rf $tmp
    mkdir -p $tmp
    ensure_offline

    local mode="full"
    [[ ! -s "$wd/augustus.hints.aa" ]] && mode="no_augustus"

    log_f "  [$idx/$n_total] [$sp] 启动 mode=$mode"
    local t0=$(date +%s)
    local rc=0

    if [[ "$mode" == "full" ]]; then
        python3 $BBC_SCRIPT \
            -m $tmp \
            -c $COMPLEASM_PATH/compleasm.py \
            -d $wd \
            -y $TSEBRA \
            -f $GETANNO \
            -g $genome \
            -t 8 \
            -p $BUSCO_LIB \
            > $log 2>&1
        rc=$?
    else
        {
          echo "# Downgraded mode: augustus.hints.aa not available for this species."
          echo "# BRAKER-internal best_by_compleasm.py requires all three gene sets;"
          echo "# here braker.aa and genemark.aa are scored with the same compleasm command."
        } > $log
        if [[ -s "$wd/GeneMark-EP/genemark.gtf" ]]; then
            python3 $GETANNO -g $genome -f $wd/GeneMark-EP/genemark.gtf -o $tmp/genemark -s True >> $log 2>&1
        fi
        python3 $COMPLEASM_PATH/compleasm.py protein -p $wd/braker.aa -l $BUSCO_LIB -t 8 -o $tmp/braker >> $log 2>&1
        rc=$?
        if [[ -s "$tmp/genemark.aa" ]]; then
            python3 $COMPLEASM_PATH/compleasm.py protein -p $tmp/genemark.aa -l $BUSCO_LIB -t 8 -o $tmp/genemark >> $log 2>&1
        fi
    fi

    local t1=$(date +%s)
    local elapsed=$((t1 - t0))
    ensure_offline

    local bM=$(getm $tmp/braker/summary.txt);   [[ -z "$bM" ]] && bM="NA"
    local aM=$(getm $tmp/augustus/summary.txt); [[ -z "$aM" ]] && aM="NA"
    local gM=$(getm $tmp/genemark/summary.txt); [[ -z "$gM" ]] && gM="NA"
    local bS=$(gets $tmp/braker/summary.txt);   [[ -z "$bS" ]] && bS="NA"
    local bD=$(getd $tmp/braker/summary.txt);   [[ -z "$bD" ]] && bD="NA"
    local bF=$(getf $tmp/braker/summary.txt);   [[ -z "$bF" ]] && bF="NA"

    if [[ "$rc" -ne 0 || "$bM" == "NA" ]]; then
        log_f "  [$idx/$n_total] [$sp] FAIL rc=$rc elapsed=${elapsed}s"
        tail -4 $log | sed 's/^/      /' | tee -a $FULL_LOG
        echo "FAIL rc=$rc" > $st
        echo -e "$sp\tFAIL_rc$rc\t$mode\t$bM\t$aM\t$gM\t$bS\t$bD\t$bF\tno\tNA\tNA\t$elapsed" >> $FULL_MANIFEST
        return 1
    fi

    local decision="unknown"
    if [[ "$mode" == "no_augustus" ]]; then
        decision="no_augustus_scored_only"
    elif grep -q "is the best one" $log; then
        decision="keep_braker_best"
    elif grep -q "there are no BUSCOs to be added" $log; then
        decision="keep_braker_no_add"
    elif grep -q "will be kept" $log; then
        decision="keep_braker"
    elif grep -qiE "merge|added to the braker" $log; then
        decision="merged"
    fi

    local final_gen="no" final_n="NA"
    if [[ -s "$wd/braker_final.gtf" ]]; then
        final_gen="yes"
        [[ -s "$wd/braker_final.aa" ]] && final_n=$(grep -c '^>' $wd/braker_final.aa)
    fi

    log_f "  [$idx/$n_total] [$sp] OK br=${bM}% au=${aM}% gm=${gM}% dec=$decision final=$final_gen ${elapsed}s"
    echo "DONE $(ts) mode=$mode braker_M=$bM decision=$decision" > $st
    echo -e "$sp\tDONE\t$mode\t$bM\t$aM\t$gM\t$bS\t$bD\t$bF\t$final_gen\t$final_n\t$decision\t$elapsed" >> $FULL_MANIFEST
    return 0
}

i=0
while read sp; do
    [[ -z "$sp" ]] && continue
    i=$((i + 1))
    run_one "$sp" "$i"
done < $SPECIES_LIST_FILE

log_f ""
log_f "===== 06c 全量完成 ====="
n_done=$(grep -l '^DONE' $STATE_DIR/*.state 2>/dev/null | wc -l)
n_fail=$(grep -l '^FAIL' $STATE_DIR/*.state 2>/dev/null | wc -l)
n_skip=$(grep -l '^SKIP' $STATE_DIR/*.state 2>/dev/null | wc -l)
n_final=$(awk -F'\t' 'NR>1 && $10=="yes"' $FULL_MANIFEST | wc -l)
log_f "  DONE:  $n_done"
log_f "  FAIL:  $n_fail"
log_f "  SKIP:  $n_skip"
log_f "  生成 braker_final.gtf: $n_final"
log_f ""
log_f "决策分布:"
awk -F'\t' 'NR>1 {c[$12]++} END {for(k in c) printf "  %-26s %d\n", k, c[k]}' $FULL_MANIFEST | tee -a $FULL_LOG
