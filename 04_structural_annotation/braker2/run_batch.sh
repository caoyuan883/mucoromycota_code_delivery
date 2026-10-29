#!/bin/bash
# ============================================================================
# Step 04 — BRAKER2 v3 Supervisor（4 并行调度，修复版）
# ============================================================================
# 修复内容:
#   ✅ count_state 精确匹配（原来 grep -c "^DONE" 会把 DONE_PARTIAL 也算成 DONE）
#   ✅ 断点续跑（73 已完成 state 自动跳过）
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

mkdir -p "$LOG_DIR" "$STATE_DIR" "$LOCK_DIR"

SUPERVISOR_LOG=$BATCH_ROOT/_supervisor.log
STOP_FLAG=$BATCH_ROOT/_STOP

log_sup() { echo "[$(ts)] [SUPV] $*" | tee -a "$SUPERVISOR_LOG"; }

if [[ ! -s "$SPECIES_LIST_FILE" ]]; then
    log_sup "ERROR: 物种列表不存在 $SPECIES_LIST_FILE"
    exit 1
fi

TOTAL=$(wc -l < "$SPECIES_LIST_FILE")
log_sup "===== BRAKER2 v3 Supervisor 启动（修复版）====="
log_sup "总物种数: $TOTAL"
log_sup "并发上限: $MAX_PARALLEL"
log_sup "线程/株:  $THREADS_PER_JOB"
log_sup "输出目录: $BATCH_ROOT"
log_sup "优雅停止: touch $STOP_FLAG"

# ---- 统计函数（精确匹配版）-----------------------------------------------
count_state() {
    local target="$1"
    local count=0
    for f in "$STATE_DIR"/*.state; do
        [[ -e "$f" ]] || continue
        local st=$(awk 'NR==1{print $1}' "$f" 2>/dev/null)
        if [[ "$target" == "FAIL" ]]; then
            [[ "$st" =~ ^FAIL ]] && count=$((count + 1))
        else
            [[ "$st" == "$target" ]] && count=$((count + 1))
        fi
    done
    echo $count
}

count_running_procs() {
    pgrep -f "run_braker2_one.sh" 2>/dev/null | wc -l
}

next_species() {
    while read -r sp; do
        [[ -z "$sp" ]] && continue
        st_file="$STATE_DIR/${sp}.state"
        if [[ ! -f "$st_file" ]]; then
            echo "$sp"
            return 0
        fi
        st=$(awk 'NR==1{print $1}' "$st_file")
        case "$st" in
            DONE|DONE_PARTIAL|RUNNING|FAIL*) continue ;;
            *) echo "$sp"; return 0 ;;
        esac
    done < "$SPECIES_LIST_FILE"
    return 1
}

tick=0
consecutive_idle=0

while true; do
    tick=$((tick + 1))
    
    if [[ -f "$STOP_FLAG" ]]; then
        log_sup "检测到 STOP flag，停止调度新任务"
        log_sup "如需重启，删除 $STOP_FLAG 后重新启动 supervisor"
        break
    fi
    
    done_ok=$(count_state "DONE")
    done_partial=$(count_state "DONE_PARTIAL")
    running=$(count_state "RUNNING")
    fail=$(count_state "FAIL")
    procs=$(count_running_procs)
    
    total_done=$((done_ok + done_partial))
    total_finished=$((total_done + fail))
    pct=$(( total_finished * 100 / TOTAL ))
    
    log_sup "tick=$tick | DONE=$done_ok DONE_PARTIAL=$done_partial RUNNING=$running FAIL=$fail | procs=$procs | ${total_finished}/${TOTAL} (${pct}%)"
    
    if [[ "$total_finished" -ge "$TOTAL" ]]; then
        log_sup "===== 全部 $TOTAL 株已完成（DONE=$total_done FAIL=$fail），supervisor 退出 ====="
        break
    fi
    
    slots_free=$((MAX_PARALLEL - procs))
    if [[ "$slots_free" -le 0 ]]; then
        consecutive_idle=0
        sleep 30
        continue
    fi
    
    launched=0
    while [[ "$launched" -lt "$slots_free" ]]; do
        sp=$(next_species || true)
        if [[ -z "$sp" ]]; then
            break
        fi
        log_sup "启动 $sp"
        nohup bash "$here/run_braker2_one.sh" "$sp" \
            > "$LOG_DIR/${sp}.sup.out" 2>&1 &
        launched=$((launched + 1))
        sleep 3
    done
    
    if [[ "$launched" -eq 0 ]]; then
        consecutive_idle=$((consecutive_idle + 1))
        if [[ "$procs" -gt 0 ]]; then
            log_sup "无新任务可启动，等已跑的 $procs 个进程完成..."
        else
            log_sup "⚠️  无新任务且无进程在跑（$total_finished/$TOTAL 完成）"
            if [[ "$consecutive_idle" -ge 5 ]]; then
                log_sup "===== 连续 5 次 idle 且无进程，supervisor 退出 ====="
                break
            fi
        fi
    else
        consecutive_idle=0
    fi
    
    sleep 30
done

log_sup ""
log_sup "===== Supervisor 结束 ====="
done_ok=$(count_state "DONE")
done_partial=$(count_state "DONE_PARTIAL")
running=$(count_state "RUNNING")
fail=$(count_state "FAIL")
log_sup "  DONE:         $done_ok"
log_sup "  DONE_PARTIAL: $done_partial"
log_sup "  RUNNING（残留）: $running"
log_sup "  FAIL:         $fail"
log_sup "  合计:         $((done_ok + done_partial + running + fail)) / $TOTAL"
log_sup ""
log_sup "详细日志: $SUPERVISOR_LOG"
