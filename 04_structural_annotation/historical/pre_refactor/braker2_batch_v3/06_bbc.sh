#!/bin/bash
# ============================================================================
# Step 06 — BUSCO/compleasm 补打分（125 株 braker.aa 统一评估）
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

BBC_ROOT=$BATCH_ROOT/_bbc
BBC_LOG_DIR=$BBC_ROOT/_logs
BBC_STATE_DIR=$BBC_ROOT/_state
BBC_SUMMARY=$BATCH_ROOT/_bbc_summary.tsv

mkdir -p $BBC_ROOT $BBC_LOG_DIR $BBC_STATE_DIR

MASTER_LOG=$BATCH_ROOT/_bbc_master.log
: > $MASTER_LOG

log_bbc() { echo "[$(ts)] [BBC] $*" | tee -a $MASTER_LOG; }

log_bbc "===== BBC 补打分启动 ====="
log_bbc "输入: $BATCH_ROOT/*/braker.aa"
log_bbc "输出: $BBC_ROOT/<species>/summary.txt"
log_bbc "汇总: $BBC_SUMMARY"
log_bbc "BUSCO 库: $BUSCO_LINEAGE_DIR"
log_bbc "并发: $MAX_PARALLEL × $THREADS_PER_JOB threads"

score_one() {
    local sp=$1
    local wd=$BATCH_ROOT/$sp
    local aa=$wd/braker.aa
    local out=$BBC_ROOT/$sp
    local state=$BBC_STATE_DIR/${sp}.state
    local stdout=$BBC_LOG_DIR/${sp}.stdout
    local stderr=$BBC_LOG_DIR/${sp}.stderr
    
    if [[ -f "$state" ]]; then
        local st=$(awk 'NR==1{print $1}' "$state")
        if [[ "$st" == "DONE" ]]; then
            return 0
        fi
    fi
    
    if [[ ! -s "$aa" ]]; then
        echo "FAIL_NO_INPUT $(ts)" > "$state"
        log_bbc "  ❌ $sp: braker.aa 不存在"
        return 1
    fi
    
    mkdir -p "$out"
    echo "RUNNING $(ts)" > "$state"
    
    local t0=$(date +%s)
    
    set +e
    python3 $COMPLEASM_PATH/compleasm.py protein \
        -p "$aa" \
        -l mucoromycota_odb12 \
        -L "$BUSCO_LINEAGE_DIR" \
        -t $THREADS_PER_JOB \
        -o "$out" \
        > "$stdout" 2> "$stderr"
    local rc=$?
    set -e
    
    local t1=$(date +%s)
    local elapsed=$((t1 - t0))
    
    local summary_file=""
    if [[ -f "$out/summary.txt" ]]; then
        summary_file="$out/summary.txt"
    elif [[ -f "$out/mucoromycota_odb12/summary.txt" ]]; then
        summary_file="$out/mucoromycota_odb12/summary.txt"
    fi
    
    if [[ -n "$summary_file" ]] && [[ -s "$summary_file" ]]; then
        echo "DONE $(ts) elapsed=${elapsed}s rc=$rc summary=$summary_file" > "$state"
        log_bbc "  ✅ $sp: elapsed=${elapsed}s"
        return 0
    else
        echo "FAIL $(ts) elapsed=${elapsed}s rc=$rc" > "$state"
        log_bbc "  ❌ $sp: rc=$rc elapsed=${elapsed}s"
        return 1
    fi
}

current_running() {
    jobs -r | wc -l
}

wait_slot() {
    while [[ $(current_running) -ge $MAX_PARALLEL ]]; do
        sleep 3
        wait -n 2>/dev/null || true
    done
}

i=0
n_total=$(wc -l < $SPECIES_LIST_FILE)
log_bbc "开始处理 $n_total 株..."

for sp in $(cat $SPECIES_LIST_FILE); do
    i=$((i + 1))
    
    state=$BBC_STATE_DIR/${sp}.state
    if [[ -f "$state" ]] && [[ "$(awk 'NR==1{print $1}' $state)" == "DONE" ]]; then
        log_bbc "[$i/$n_total] $sp: 已完成，跳过"
        continue
    fi
    
    wait_slot
    log_bbc "[$i/$n_total] 启动 $sp"
    score_one "$sp" &
done

wait
log_bbc "===== 全部完成 ====="

log_bbc "生成汇总表: $BBC_SUMMARY"

echo -e "species\tstate\tC_single\tC_dup\tF\tM\ttotal\tproteins\tsource_summary" > $BBC_SUMMARY

for sp in $(cat $SPECIES_LIST_FILE); do
    state=$BBC_STATE_DIR/${sp}.state
    aa=$BATCH_ROOT/$sp/braker.aa
    
    n_prot="NA"
    [[ -s "$aa" ]] && n_prot=$(grep -c '^>' "$aa")
    
    st="MISSING"
    [[ -f "$state" ]] && st=$(awk 'NR==1{print $1}' "$state")
    
    summary=""
    for f in \
        "$BBC_ROOT/$sp/summary.txt" \
        "$BBC_ROOT/$sp/mucoromycota_odb12/summary.txt"; do
        [[ -s "$f" ]] && summary=$f && break
    done
    
    if [[ -z "$summary" ]]; then
        echo -e "$sp\t$st\tNA\tNA\tNA\tNA\tNA\t$n_prot\tNONE" >> $BBC_SUMMARY
        continue
    fi
    
    c_s=$(awk -F'[, %]+' '/^S:/{print $2}' $summary | head -1)
    c_d=$(awk -F'[, %]+' '/^D:/{print $2}' $summary | head -1)
    c_f=$(awk -F'[, %]+' '/^F:/{print $2}' $summary | head -1)
    c_m=$(awk -F'[, %]+' '/^M:/{print $2}' $summary | head -1)
    c_n=$(awk -F'[, ]+' '/^N:/{print $2}' $summary | head -1)
    
    echo -e "$sp\t$st\t${c_s:-NA}\t${c_d:-NA}\t${c_f:-NA}\t${c_m:-NA}\t${c_n:-NA}\t$n_prot\t$summary" >> $BBC_SUMMARY
done

log_bbc "===== 汇总表已生成: $BBC_SUMMARY ====="

n_done=$(grep -l "^DONE" $BBC_STATE_DIR/*.state 2>/dev/null | wc -l)
n_fail=$(grep -l "^FAIL" $BBC_STATE_DIR/*.state 2>/dev/null | wc -l)
log_bbc ""
log_bbc "===== BBC 最终统计 ====="
log_bbc "  ✅ DONE: $n_done / $n_total"
log_bbc "  ❌ FAIL: $n_fail"
log_bbc "  汇总表: $BBC_SUMMARY"

if [[ "$n_done" -eq "$n_total" ]]; then
    log_bbc "  🎉 全部成功"
else
    log_bbc "  ⚠️  有失败，检查 $BBC_LOG_DIR/<species>.stderr"
fi
