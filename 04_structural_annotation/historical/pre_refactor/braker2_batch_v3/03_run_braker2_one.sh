#!/bin/bash
# ============================================================================
# Step 03 — 单株 BRAKER2 v3 启动
# 用法: bash 03_run_braker2_one.sh <species_name>
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

SP=${1:-}
if [[ -z "$SP" ]]; then
    echo "用法: bash 03_run_braker2_one.sh <species_name>" >&2
    exit 2
fi

mkdir -p "$BATCH_ROOT" "$LOG_DIR" "$STATE_DIR" "$LOCK_DIR"

FA=$GENOME_DIR/${SP}_genomic_softmasked.clean.fa
WORKDIR=$BATCH_ROOT/$SP
STATE=$STATE_DIR/${SP}.state
LOCK=$LOCK_DIR/${SP}.lock
STDOUT=$LOG_DIR/${SP}.stdout
STDERR=$LOG_DIR/${SP}.stderr

if [[ ! -s "$FA" ]]; then
    echo "[$SP] ❌ clean FASTA 不存在: $FA" >&2
    echo "FAIL_NO_INPUT $(ts)" > "$STATE"
    exit 3
fi
if [[ ! -s "$PROT_DB" ]]; then
    echo "[$SP] ❌ Protein DB 不存在: $PROT_DB" >&2
    echo "FAIL_NO_INPUT $(ts)" > "$STATE"
    exit 3
fi

exec 200>"$LOCK"
if ! flock -n 200; then
    echo "[$SP] ⚠️  另一个进程已持有锁，本次不启动" >&2
    exit 0
fi

if [[ -f "$STATE" ]]; then
    st=$(awk '{print $1}' "$STATE")
    case "$st" in
        DONE|DONE_PARTIAL)
            echo "[$SP] ✔  已完成 (state=$st)，跳过" >&2
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
echo "RUNNING $(ts) pid=$$" > "$STATE"

log "[$SP] 🚀 启动 BRAKER2 v3 (protein-only, --fungus + compleasm 增强)"
log "[$SP]    FASTA:      $FA"
log "[$SP]    Protein DB: $PROT_DB"
log "[$SP]    WD:         $WORKDIR"
log "[$SP]    GeneMark:   $GENEMARK_PATH"
log "[$SP]    TSEBRA:     $TSEBRA_PATH"
log "[$SP]    v3 增强:    --fungus + --busco_lineage + compleasm"

t0=$(date +%s)

set +e
perl "$BRAKER_SCRIPT" \
    --genome="$FA" \
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
    --PROTHINT_PATH="$PROTHINT_PATH" \
    --AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" \
    --AUGUSTUS_BIN_PATH="$AUGUSTUS_BIN_PATH" \
    --AUGUSTUS_SCRIPTS_PATH="$AUGUSTUS_SCRIPTS_PATH" \
    --BAMTOOLS_PATH="$BAMTOOLS_PATH" \
    --DIAMOND_PATH="$DIAMOND_PATH" \
    --TSEBRA_PATH="$TSEBRA_PATH" \
    --skip_fixing_broken_genes \
    --verbosity=4 \
    >>"$STDOUT" 2>>"$STDERR"
rc=$?
set -e

t1=$(date +%s)
elapsed=$((t1 - t0))
hh=$((elapsed / 3600)); mm=$(( (elapsed % 3600) / 60 )); ss=$((elapsed % 60))
elapsed_str=$(printf "%02d:%02d:%02d" $hh $mm $ss)

MIN_PROTEINS=5000
if [[ -s "$WORKDIR/braker.aa" ]]; then
    proteins=$(grep -c '^>' "$WORKDIR/braker.aa" || echo 0)
    if [[ "$proteins" -gt "$MIN_PROTEINS" ]]; then
        if [[ "$rc" -eq 0 ]]; then
            echo "DONE $(ts) elapsed=$elapsed_str proteins=$proteins rc=0" > "$STATE"
            log "[$SP] ✅ 完成，耗时 $elapsed_str，蛋白 $proteins 条"
        else
            echo "DONE_PARTIAL $(ts) elapsed=$elapsed_str proteins=$proteins rc=$rc" > "$STATE"
            log "[$SP] ⚠️  主结果齐全但收尾 rc=$rc，视为 DONE (蛋白 $proteins 条)"
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
