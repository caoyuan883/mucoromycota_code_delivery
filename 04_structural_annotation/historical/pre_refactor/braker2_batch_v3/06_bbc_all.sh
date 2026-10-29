#!/bin/bash
# ============================================================================
# Step 06 — BBC 全量打分（跟 BRAKER3 内部方法学一致）
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

BUSCO_LIB=/media/desk16/tly8529/mucoromycota_odb12

BBC_ROOT=$BATCH_ROOT/_bbc_all
BBC_LOG_DIR=$BBC_ROOT/_logs
BBC_STATE_DIR=$BBC_ROOT/_state
BBC_SUMMARY=$BATCH_ROOT/_bbc_summary.tsv

mkdir -p $BBC_ROOT $BBC_LOG_DIR $BBC_STATE_DIR

MASTER_LOG=$BATCH_ROOT/_bbc_master.log
: > $MASTER_LOG

log_bbc() { echo "[$(ts)] [BBC] $*" | tee -a $MASTER_LOG; }

log_bbc "===== BBC 全量打分启动（BRAKER3 方法学一致）====="
log_bbc "BUSCO 库: $BUSCO_LIB"
log_bbc "输出: $BBC_ROOT"
log_bbc "并发: $MAX_PARALLEL × $THREADS_PER_JOB threads"

score_one_set() {
    local sp=$1
    local sname=$2
    local aa=$3
    local out=$BBC_ROOT/$sp/$sname
    local state=$BBC_STATE_DIR/${sp}_${sname}.state
    local stdout=$BBC_LOG_DIR/${sp}_${sname}.stdout
    local stderr=$BBC_LOG_DIR/${sp}_${sname}.stderr
    
    if [[ -f "$state" ]]; then
        local st=$(awk 'NR==1{print $1}' "$state")
        if [[ "$st" == "DONE" ]]; then
            return 0
        fi
    fi
    
    if [[ ! -s "$aa" ]]; then
        echo "SKIP_NO_INPUT $(ts)" > "$state"
        return 0
    fi
    
    mkdir -p "$out"
    echo "RUNNING $(ts)" > "$state"
    
    local t0=$(date +%s)
    
    set +e
    python3 $COMPLEASM_PATH/compleasm.py protein \
        -p "$aa" \
        -l $BUSCO_LIB \
        -t $THREADS_PER_JOB \
        -o "$out" \
        > "$stdout" 2> "$stderr"
    local rc=$?
    set -e
    
    local t1=$(date +%s)
    local elapsed=$((t1 - t0))
    
    local summary_file=""
    for f in "$out/summary.txt" "$out/mucoromycota_odb12/summary.txt"; do
        [[ -s "$f" ]] && summary_file=$f && break
    done
    
    if [[ -n "$summary_file" ]]; then
        echo "DONE $(ts) elapsed=${elapsed}s rc=$rc summary=$summary_file" > "$state"
        return 0
    else
        echo "FAIL $(ts) elapsed=${elapsed}s rc=$rc" > "$state"
        return 1
    fi
}

wait_slot() {
    while [[ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]]; do
        sleep 3
        wait -n 2>/dev/null || true
    done
}

log_bbc "扫描 3 集合任务队列..."

TASK_QUEUE=$BBC_ROOT/_task_queue.txt
: > $TASK_QUEUE

n_braker=0; n_aug=0; n_gm=0
for sp in $(cat $SPECIES_LIST_FILE); do
    [[ -s "$BATCH_ROOT/$sp/braker.aa" ]] && { echo "$sp|braker|$BATCH_ROOT/$sp/braker.aa" >> $TASK_QUEUE; n_braker=$((n_braker+1)); }
    [[ -s "$BATCH_ROOT/$sp/augustus.hints.aa" ]] && { echo "$sp|augustus|$BATCH_ROOT/$sp/augustus.hints.aa" >> $TASK_QUEUE; n_aug=$((n_aug+1)); }
    [[ -s "$BATCH_ROOT/$sp/bbc/genemark.aa" ]] && { echo "$sp|genemark|$BATCH_ROOT/$sp/bbc/genemark.aa" >> $TASK_QUEUE; n_gm=$((n_gm+1)); }
done

n_tasks=$(wc -l < $TASK_QUEUE)
log_bbc "  braker:   $n_braker 株"
log_bbc "  augustus: $n_aug 株"
log_bbc "  genemark: $n_gm 株"
log_bbc "  总任务:   $n_tasks"

i=0
while IFS='|' read -r sp sname aa; do
    i=$((i + 1))
    
    state=$BBC_STATE_DIR/${sp}_${sname}.state
    if [[ -f "$state" ]] && [[ "$(awk 'NR==1{print $1}' $state)" == "DONE" ]]; then
        continue
    fi
    
    wait_slot
    log_bbc "[$i/$n_tasks] $sp / $sname"
    score_one_set "$sp" "$sname" "$aa" &
done < $TASK_QUEUE

wait
log_bbc "===== 全部完成，开始生成汇总 ====="

echo -e "species\tset\tstate\tC_single\tC_dup\tF\tM\ttotal_N\tproteins\tsource_summary" > $BBC_SUMMARY

for sp in $(cat $SPECIES_LIST_FILE); do
    for sname in braker augustus genemark; do
        state=$BBC_STATE_DIR/${sp}_${sname}.state
        
        case "$sname" in
            braker)   aa=$BATCH_ROOT/$sp/braker.aa ;;
            augustus) aa=$BATCH_ROOT/$sp/augustus.hints.aa ;;
            genemark) aa=$BATCH_ROOT/$sp/bbc/genemark.aa ;;
        esac
        n_prot="NA"
        [[ -s "$aa" ]] && n_prot=$(grep -c '^>' "$aa")
        
        st="MISSING"
        [[ -f "$state" ]] && st=$(awk 'NR==1{print $1}' "$state")
        
        summary=""
        for f in \
            "$BBC_ROOT/$sp/$sname/summary.txt" \
            "$BBC_ROOT/$sp/$sname/mucoromycota_odb12/summary.txt"; do
            [[ -s "$f" ]] && summary=$f && break
        done
        
        if [[ -z "$summary" ]]; then
            echo -e "$sp\t$sname\t$st\tNA\tNA\tNA\tNA\tNA\t$n_prot\tNONE" >> $BBC_SUMMARY
            continue
        fi
        
        c_s=$(awk '/^S:/{gsub(/%|,/,""); print $2}' $summary | head -1)
        c_d=$(awk '/^D:/{gsub(/%|,/,""); print $2}' $summary | head -1)
        c_f=$(awk '/^F:/{gsub(/%|,/,""); print $2}' $summary | head -1)
        c_m=$(awk '/^M:/{gsub(/%|,/,""); print $2}' $summary | head -1)
        c_n=$(awk '/^N:/{print $2}' $summary | head -1)
        
        echo -e "$sp\t$sname\t$st\t${c_s:-NA}\t${c_d:-NA}\t${c_f:-NA}\t${c_m:-NA}\t${c_n:-NA}\t$n_prot\t$summary" >> $BBC_SUMMARY
    done
done

log_bbc "===== 汇总表: $BBC_SUMMARY ====="

n_done=$(grep -l "^DONE" $BBC_STATE_DIR/*.state 2>/dev/null | wc -l)
n_fail=$(grep -l "^FAIL" $BBC_STATE_DIR/*.state 2>/dev/null | wc -l)
n_missing_input=$(grep -l "^SKIP_NO_INPUT" $BBC_STATE_DIR/*.state 2>/dev/null | wc -l)

log_bbc ""
log_bbc "===== BBC 最终统计 ====="
log_bbc "  ✅ DONE:            $n_done"
log_bbc "  ❌ FAIL:            $n_fail"
log_bbc "  ⏭️  SKIP_NO_INPUT:  $n_missing_input"
log_bbc "  预期: 351 (125 braker + 101 augustus + 125 genemark)"
log_bbc "  汇总表: $BBC_SUMMARY"
log_bbc ""
log_bbc "查看汇总（前 10 行）:"
head -10 $BBC_SUMMARY | column -t -s $'\t' | tee -a $MASTER_LOG
