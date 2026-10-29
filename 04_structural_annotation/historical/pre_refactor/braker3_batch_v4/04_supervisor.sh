#!/bin/bash
# ============================================================================
# BRAKER3 v3 调度器 — 4 并行动态调度
# 跑完一株立刻启动下一株（不等一整批）
#
# 使用：
#   nohup bash 04_supervisor.sh > _logs/SUPV.out 2>&1 &
#
# 停止：
#   pkill -f 04_supervisor.sh
#   （BRAKER 子进程会继续跑完，不受影响）
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

if [[ ${#SPECIES_BRAKER3[@]} -eq 0 ]]; then
  echo "❌ 物种表为空，请先跑 00_prep_species_list.sh" >&2
  exit 1
fi

mkdir -p "$BATCH_ROOT" "$LOG_DIR" "$STATE_DIR" "$LOCK_DIR"

SUPV_LOG=$LOG_DIR/SUPV_$(date +%Y%m%d_%H%M%S).log

log_supv() { echo "[$(ts)] $*" | tee -a "$SUPV_LOG"; }

log_supv "════════════════════════════════════════════════"
log_supv "BRAKER3 v3 supervisor 启动"
log_supv "并发数: $MAX_PARALLEL   线程/株: $THREADS_PER_JOB"
log_supv "总物种数: ${#SPECIES_BRAKER3[@]}"
log_supv "输出根: $BATCH_ROOT"
log_supv "v3 增强: --fungus + --busco_lineage + compleasm"
log_supv "════════════════════════════════════════════════"

# ---- 判断一株是否 running -------------------------------------------------
is_running() {
  local sp=$1
  local lock=$LOCK_DIR/${sp}.lock
  [[ -f "$lock" ]] || return 1
  fuser "$lock" >/dev/null 2>&1
}

# ---- 判断一株是否需要跑 ---------------------------------------------------
needs_run() {
  local sp=$1
  local state=$STATE_DIR/${sp}.state
  if is_running "$sp"; then return 1; fi
  if [[ -f "$state" ]]; then
    local st=$(awk '{print $1}' "$state")
    case "$st" in
      DONE)          return 1 ;;
      DONE_PARTIAL)  return 1 ;;
      RUNNING)       return 1 ;;
      FAIL*)         return 1 ;;
    esac
  fi
  return 0
}

count_running() {
  local n=0
  for sp in "${SPECIES_BRAKER3[@]}"; do
    if is_running "$sp"; then n=$((n+1)); fi
  done
  echo "$n"
}

count_done() {
  local n=0
  for sp in "${SPECIES_BRAKER3[@]}"; do
    local state=$STATE_DIR/${sp}.state
    if [[ -f "$state" ]] && awk '{print $1}' "$state" | grep -q '^DONE$'; then
      n=$((n+1))
    fi
  done
  echo "$n"
}

count_fail() {
  local n=0
  for sp in "${SPECIES_BRAKER3[@]}"; do
    local state=$STATE_DIR/${sp}.state
    if [[ -f "$state" ]] && awk '{print $1}' "$state" | grep -q '^FAIL'; then
      n=$((n+1))
    fi
  done
  echo "$n"
}

# ---- 快速失败探测（两道保险）-----------------------------------------------
# 保险 A：启动一株后强制等 60 秒再判断 running（避免 BRAKER 未来得及创建 lock 就又启新株）
# 保险 B：检查最新 3 株是否都在 <5 分钟内 FAIL（环境破了 → 立即停）
check_rapid_fail() {
  # 列出最近 3 株有 state 的，看是否全部 FAIL 且 elapsed < 300秒
  local recent_states
  recent_states=$(ls -t "$STATE_DIR"/*.state 2>/dev/null | head -3)
  local n_files
  n_files=$(echo "$recent_states" | grep -c "." || echo 0)
  if [[ "$n_files" -lt 3 ]]; then return 0; fi

  local all_fast_fail=1
  for f in $recent_states; do
    local line
    line=$(cat "$f")
    local status
    status=$(echo "$line" | awk '{print $1}')
    # 不是 FAIL 就不算
    if [[ "$status" != FAIL* ]]; then all_fast_fail=0; break; fi
    # 提 elapsed=HH:MM:SS
    local elapsed_str
    elapsed_str=$(echo "$line" | grep -oE 'elapsed=[0-9:]+' | cut -d= -f2)
    if [[ -z "$elapsed_str" ]]; then all_fast_fail=0; break; fi
    local hh mm ss total
    IFS=: read -r hh mm ss <<< "$elapsed_str"
    total=$((10#$hh*3600 + 10#$mm*60 + 10#$ss))
    if [[ "$total" -ge 300 ]]; then all_fast_fail=0; break; fi
  done

  if [[ "$all_fast_fail" -eq 1 ]]; then
    log_supv "❌ 紧急止损：最新 3 株都在 5 分钟内 FAIL，环境可能有问题！"
    log_supv "    最近 3 株状态:"
    for f in $recent_states; do
      log_supv "    - $(basename $f .state): $(cat $f)"
    done
    log_supv "    supervisor 主动退出。请检查日志后手已重启。"
    return 1
  fi
  return 0
}

# ---- 主循环 ---------------------------------------------------------------
TOTAL=${#SPECIES_BRAKER3[@]}
while true; do
  running=$(count_running)
  done_n=$(count_done)
  fail_n=$(count_fail)
  pending=$(( TOTAL - running - done_n - fail_n ))

  # 全部结束？
  if [[ "$done_n" -eq "$TOTAL" ]]; then
    log_supv "🎉 全部 $TOTAL 株 DONE，调度器退出"
    break
  fi
  if [[ "$running" -eq 0 && "$pending" -eq 0 ]]; then
    log_supv "⚠️  没有 running 也没有 pending  (DONE=$done_n  FAIL=$fail_n)"
    log_supv "    如需重跑 FAIL 株，请先 \\rm 对应 state + 移走 workdir，然后重启 supervisor"
    break
  fi

  # 保险 B：在启动新株前先检查快速失败
  if ! check_rapid_fail; then
    break
  fi

  # 有空闲位置就启动新的
  while [[ "$running" -lt "$MAX_PARALLEL" ]]; do
    launched=0
    for sp in "${SPECIES_BRAKER3[@]}"; do
      if needs_run "$sp"; then
        log_supv "[LAUNCH] $sp  (当前 running=$running  DONE=$done_n)"
        nohup bash "$here/03_run_braker3_one.sh" "$sp" \
          >/dev/null 2>>"$LOG_DIR/${sp}.launch.err" &
        launched=1
        # 保险 A：启动后强制等 60 秒（避免 BRAKER 未创建 lock 就又启新株）
        log_supv "    启动后等待 60 秒确认 workdir/lock 建立……"
        sleep 60
        # 重新检查快速失败（如果刚启动的那株已经 FAIL了）
        if ! check_rapid_fail; then
          exit 1
        fi
        running=$(count_running)
        break
      fi
    done
    if [[ "$launched" -eq 0 ]]; then break; fi
    if [[ "$running" -ge "$MAX_PARALLEL" ]]; then break; fi
  done

  # 简报
  running=$(count_running)
  done_n=$(count_done)
  fail_n=$(count_fail)
  pending=$(( TOTAL - running - done_n - fail_n ))
  log_supv "STATUS  running=$running  DONE=$done_n  FAIL=$fail_n  pending=$pending"

  sleep 30
done

log_supv "════════════════════════════════════════════════"
log_supv "supervisor 结束"
log_supv "最终: DONE=$(count_done)  FAIL=$(count_fail)  RUNNING=$(count_running)"
log_supv "════════════════════════════════════════════════"
