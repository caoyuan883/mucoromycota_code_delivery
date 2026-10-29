#!/bin/bash
# ============================================================================
# Step 14 — 新增 2 株 BRAKER2（方法学与 125 株完全一致）
#   1) header 清洗 → 03_softmasked_genomes_v3/<sp>_genomic_softmasked.clean.fa
#   2) 调用原班 03_run_braker2_one.sh，参数一字未改
#   3) 8 线程 × 2 并行
# 安全: 绝不删除任何 workdir；已存在则拒绝启动
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

NEWDIR=$DATA_ROOT/03_softmasked_genomes_260807
SPECIES="Cryptocalciella_humilis Rhizopus_sp._TYR1"
LOG=$BATCH_ROOT/_addon_260807.log

mkdir -p "$LOG_DIR" "$STATE_DIR" "$LOCK_DIR"

echo "[$(ts)] ===== 新增 2 株 BRAKER2 =====" | tee -a $LOG

# ---------------- Part 1: header 清洗 ----------------
echo "[$(ts)] --- header 清洗 ---" | tee -a $LOG
for sp in $SPECIES; do
    raw=$NEWDIR/${sp}_genomic_softmasked.fa
    clean=$GENOME_DIR/${sp}_genomic_softmasked.clean.fa

    if [ ! -s "$raw" ]; then
        echo "  [$sp] 原始文件缺失: $raw" | tee -a $LOG; exit 3
    fi
    if [ -s "$clean" ]; then
        echo "  [$sp] clean.fa 已存在，跳过清洗" | tee -a $LOG; continue
    fi

    awk '/^>/{print $1; next} {print}' "$raw" > "$clean.tmp"
    rl=$(wc -l < "$raw"); tl=$(wc -l < "$clean.tmp")
    if [ "$rl" != "$tl" ]; then
        echo "  [$sp] 行数不符 raw=$rl tmp=$tl，放弃" | tee -a $LOG
        \rm -f "$clean.tmp"; exit 4
    fi
    if head -1 "$clean.tmp" | grep -q ' '; then
        echo "  [$sp] header 仍含空格，放弃" | tee -a $LOG
        \rm -f "$clean.tmp"; exit 5
    fi
    mv "$clean.tmp" "$clean"
    echo "  [$sp] OK  $(du -h $clean | cut -f1)  contigs=$(grep -c '^>' $clean)  首行=$(head -1 $clean)" | tee -a $LOG
done

# ---------------- Part 2: 预检 ----------------
echo "" | tee -a $LOG
echo "[$(ts)] --- 预检 ---" | tee -a $LOG
[ -s "$PROT_DB" ] && echo "  蛋白库 OK: $PROT_DB ($(du -h $PROT_DB | cut -f1))" | tee -a $LOG || { echo "  蛋白库缺失"; exit 6; }
[ -s "$BRAKER_SCRIPT" ] && echo "  braker.pl OK" | tee -a $LOG || { echo "  braker.pl 缺失"; exit 6; }
[ -d "$BUSCO_LINEAGE_DIR" ] && echo "  BUSCO 库 OK: $BUSCO_LINEAGE_DIR" | tee -a $LOG || { echo "  BUSCO 库缺失"; exit 6; }
echo "  可用空间: $(df -h /media/desk16 | awk 'NR==2{print $4}')" | tee -a $LOG

for sp in $SPECIES; do
    [ -e "$BATCH_ROOT/$sp" ] && { echo "  [$sp] workdir 已存在，拒绝启动（本脚本不删任何目录）" | tee -a $LOG; exit 7; }
    [ -e "$STATE_DIR/${sp}.state" ] && { echo "  [$sp] state 已存在: $(cat $STATE_DIR/${sp}.state)" | tee -a $LOG; exit 7; }
    [ -e "$AUGUSTUS_CONFIG_PATH/species/${sp}_${SPECIES_SUFFIX}" ] && \
        echo "  [$sp] 注意: augustus species 目录已存在，BRAKER 会复用" | tee -a $LOG
done

# ---------------- Part 3: 启动，8 线程 × 2 并行 ----------------
echo "" | tee -a $LOG
echo "[$(ts)] --- 启动 BRAKER2（8 线程 × 2 并行）---" | tee -a $LOG
for sp in $SPECIES; do
    nohup bash "$here/03_run_braker2_one.sh" "$sp" >> "$LOG_DIR/${sp}.launch.out" 2>&1 &
    echo "  [$sp] 已启动 pid=$!" | tee -a $LOG
    sleep 5
done

echo "" | tee -a $LOG
echo "[$(ts)] 两株均已在后台运行。监控：" | tee -a $LOG
echo "  cat $STATE_DIR/Cryptocalciella_humilis.state" | tee -a $LOG
echo "  cat $STATE_DIR/Rhizopus_sp._TYR1.state" | tee -a $LOG
echo "  tail -f $LOG_DIR/Cryptocalciella_humilis.stderr" | tee -a $LOG
wait
echo "[$(ts)] ===== 两株全部结束 =====" | tee -a $LOG
for sp in $SPECIES; do
    echo "  [$sp] $(cat $STATE_DIR/${sp}.state 2>/dev/null)" | tee -a $LOG
done
