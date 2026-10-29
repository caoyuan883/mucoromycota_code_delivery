#!/bin/bash
# ============================================================================
# 单株 BRAKER3 启动 v3 — Scientific Data 正式版
# 用法: bash run_braker3_one.sh <species_name>
#
# vs v2 关键差异：
#   ✅ 使用 _genomic_softmasked.clean.fa（Step 0 清洗后）
#   ✅ 新增 --fungus
#   ✅ 新增 --busco_lineage
#   ✅ 新增 --COMPLEASM_PATH
#   ✅ 物种名后缀 _braker3_v3
#
# 保持 v2 已验证的防灾机制：
#   ❌ 绝不 rm -rf workdir
#   ✅ flock 互斥锁
#   ✅ 状态严格化
#   ✅ 双写 state
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

SP=${1:-}
if [[ -z "$SP" ]]; then
  echo "用法: bash run_braker3_one.sh <species_name>" >&2
  exit 2
fi

mkdir -p "$BATCH_ROOT" "$LOG_DIR" "$STATE_DIR" "$LOCK_DIR"

FA=$GENOME_DIR/${SP}_genomic_softmasked.clean.fa
BAM=$BAM_DIR/${SP}.merged.bam
WORKDIR=$BATCH_ROOT/$SP
STATE=$STATE_DIR/${SP}.state
LOCK=$LOCK_DIR/${SP}.lock
STDOUT=$LOG_DIR/${SP}.stdout
STDERR=$LOG_DIR/${SP}.stderr

# ---- 输入检查 --------------------------------------------------------------
if [[ ! -s "$FA" ]]; then
  echo "[$SP] ❌ clean FASTA 不存在: $FA" >&2
  echo "FAIL_NO_INPUT $(ts)" > "$STATE"
  exit 3
fi
if [[ ! -s "$BAM" ]]; then
  echo "[$SP] ❌ BAM 不存在: $BAM" >&2
  echo "FAIL_NO_INPUT $(ts)" > "$STATE"
  exit 3
fi

# ---- 获取锁 ----------------------------------------------------------------
exec 200>"$LOCK"
if ! flock -n 200; then
  echo "[$SP] ⚠️  另一个进程已持有锁，本次不启动" >&2
  exit 0
fi

# ---- 状态检查（严格）-----------------------------------------------------
if [[ -f "$STATE" ]]; then
  st=$(awk '{print $1}' "$STATE")
  case "$st" in
    DONE)
      echo "[$SP] ✔  已完成，跳过" >&2
      exit 0
      ;;
    RUNNING)
      echo "[$SP] ⚠️  state=RUNNING，可能是残留脏状态；本进程不重启，请人工检查" >&2
      exit 4
      ;;
    FAIL*)
      echo "[$SP] ❌ state=$st，需要人工先重置状态才能重跑：" >&2
      echo "        \\rm $STATE" >&2
      echo "        mv $WORKDIR ${WORKDIR}.bak_\$(date +%s)  # 若存在" >&2
      exit 5
      ;;
    *)
      echo "[$SP] state 内容异常: $st，拒绝启动" >&2
      exit 6
      ;;
  esac
fi

# ---- Workdir 检查（绝对不 rm）---------------------------------------------
if [[ -d "$WORKDIR" ]]; then
  count=$(ls -A "$WORKDIR" 2>/dev/null | wc -l)
  if [[ "$count" -gt 0 ]]; then
    echo "[$SP] ❌ workdir 已存在且非空 ($count 项): $WORKDIR" >&2
    echo "        本脚本绝不删除任何 workdir。请人工确认后手动处理：" >&2
    echo "          mv $WORKDIR ${WORKDIR}.bak_\$(date +%s)" >&2
    echo "          \\rm $STATE 2>/dev/null" >&2
    exit 7
  fi
fi

mkdir -p "$WORKDIR"

# ---- 写 state = RUNNING ----------------------------------------------------
echo "RUNNING $(ts) pid=$$" > "$STATE"

# ---- 启动 BRAKER3 v3 -------------------------------------------------------
log "[$SP] 🚀 启动 BRAKER3 v3"
log "[$SP]    FASTA: $FA"
log "[$SP]    BAM:   $BAM"
log "[$SP]    WD:    $WORKDIR"
log "[$SP]    v3 增强: --fungus + --busco_lineage + compleasm"

t0=$(date +%s)

set +e
perl "$BRAKER_SCRIPT" \
  --genome="$FA" \
  --bam="$BAM" \
  --prot_seq="$PROT_DB" \
  --softmasking \
  --species=${SP}_${SPECIES_SUFFIX} \
  --workingdir="$WORKDIR" \
  --threads=$THREADS_PER_JOB \
  --gff3 \
  --fungus \
  --busco_lineage="$BUSCO_LINEAGE_DIR" \
  --COMPLEASM_PATH="$COMPLEASM_PATH" \
  --GENEMARK_PATH="$GENEMARK_PATH" \
  --BAMTOOLS_PATH="$BAMTOOLS_PATH" \
  --DIAMOND_PATH="$DIAMOND_PATH" \
  --PROTHINT_PATH="$PROTHINT_PATH" \
  --TSEBRA_PATH="$TSEBRA_PATH" \
  --skip_fixing_broken_genes \
  >>"$STDOUT" 2>>"$STDERR"
rc=$?
set -e

t1=$(date +%s)
elapsed=$((t1 - t0))
hh=$((elapsed / 3600)); mm=$(( (elapsed % 3600) / 60 )); ss=$((elapsed % 60))
elapsed_str=$(printf "%02d:%02d:%02d" $hh $mm $ss)

# ---- 结果判定 --------------------------------------------------------------
# v3 策略：主产物 braker.aa 存在且非空即视为 DONE
# 即使 rc != 0（例如 best_by_compleasm 收尾步骤失败），只要主结果齐全就接受
# 蛋白数需 > 5000（防止空文件误判）
MIN_PROTEINS=5000
if [[ -s "$WORKDIR/braker.aa" ]]; then
  proteins=$(grep -c '^>' "$WORKDIR/braker.aa" || echo 0)
  if [[ "$proteins" -gt "$MIN_PROTEINS" ]]; then
    if [[ "$rc" -eq 0 ]]; then
      echo "DONE $(ts) elapsed=$elapsed_str proteins=$proteins rc=0" > "$STATE"
      log "[$SP] ✅ 完成，耗时 $elapsed_str，蛋白 $proteins 条"
    else
      # rc != 0 但主产物齐全 —— best_by_compleasm 等收尾步骤失败
      echo "DONE_PARTIAL $(ts) elapsed=$elapsed_str proteins=$proteins rc=$rc" > "$STATE"
      log "[$SP] ⚠️  主结果齐全但收尾步骤 rc=$rc，仍视为 DONE (蛋白 $proteins 条)"
    fi
    exit 0
  else
    echo "FAIL_LOW_PROTEIN $(ts) elapsed=$elapsed_str proteins=$proteins rc=$rc" > "$STATE"
    log "[$SP] ❌ 蛋白数 $proteins 过低 (< $MIN_PROTEINS)  查看 $STDERR"
    exit 1
  fi
else
  echo "FAIL $(ts) elapsed=$elapsed_str rc=$rc" > "$STATE"
  log "[$SP] ❌ 失败 rc=$rc  耗时 $elapsed_str  查看 $STDERR"
  exit "$rc"
fi
