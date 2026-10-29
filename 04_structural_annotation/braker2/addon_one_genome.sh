#!/bin/bash
# ============================================================================
# Step 17 — 补跑单株 BRAKER2（方法学与已完成的 127 株完全一致）
# 用法: bash addon_one_genome.sh <species_name>
#   1) 自动定位软屏蔽基因组
#   2) header 清洗 → 03_softmasked_genomes_v3/<sp>_genomic_softmasked.clean.fa
#   3) 追加进物种清单（查重）
#   4) 调用原班 run_braker2_one.sh，参数一字未改
# 安全: 绝不删除任何 workdir；已存在则拒绝启动
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

SP=${1:-}
[[ -z "$SP" ]] && { echo "用法: bash addon_one_genome.sh <species_name>"; exit 2; }

LOG=$BATCH_ROOT/_addon_${SP}.log
mkdir -p "$LOG_DIR" "$STATE_DIR" "$LOCK_DIR"

echo "[$(ts)] ===== 补跑 $SP =====" | tee -a $LOG

# ---------------- 1. 定位原始软屏蔽基因组 ----------------
RAW=""
for cand in \
    "$RAW_GENOME_DIR/${SP}_genomic_softmasked.fa" \
    "$DATA_ROOT/03_softmasked_genomes/${SP}_genomic_softmasked.fa" \
    "$DATA_ROOT/03_softmasked_genomes_260807/${SP}_genomic_softmasked.fa" \
    "$DATA_ROOT/03_softmasked_genomes_v3/${SP}_genomic_softmasked.fa"
do
    [ -s "$cand" ] && { RAW=$cand; break; }
done
if [ -z "$RAW" ]; then
    echo "  未找到软屏蔽基因组，尝试全盘搜索 ..." | tee -a $LOG
    RAW=$(find $DATA_ROOT -maxdepth 2 -name "${SP}_genomic_softmasked.fa" -size +1M 2>/dev/null | head -1)
fi
[ -z "$RAW" ] && { echo "  找不到 ${SP}_genomic_softmasked.fa，请人工指定路径后再跑" | tee -a $LOG; exit 3; }
echo "  原始基因组: $RAW  ($(du -h $RAW | cut -f1))" | tee -a $LOG

# ---------------- 2. 冲突检查（绝不覆盖任何东西）----------------
if [ -e "$BATCH_ROOT/$SP" ]; then
    echo "  BRAKER2 工作目录已存在，拒绝启动: $BATCH_ROOT/$SP" | tee -a $LOG
    echo "  如确需重跑，请人工执行: mv $BATCH_ROOT/$SP ${BATCH_ROOT}/${SP}.bak_\$(date +%s)" | tee -a $LOG
    exit 7
fi
if [ -e "$STATE_DIR/${SP}.state" ]; then
    echo "  state 已存在: $(cat $STATE_DIR/${SP}.state)" | tee -a $LOG
    echo "  如确需重跑，请人工执行: \\rm $STATE_DIR/${SP}.state" | tee -a $LOG
    exit 7
fi
if [ -e "$AUGUSTUS_CONFIG_PATH/species/${SP}_${SPECIES_SUFFIX}" ]; then
    echo "  注意: augustus species 目录已存在，BRAKER 会复用：${SP}_${SPECIES_SUFFIX}" | tee -a $LOG
fi

# ---------------- 3. header 清洗 ----------------
CLEAN=$GENOME_DIR/${SP}_genomic_softmasked.clean.fa
if [ -s "$CLEAN" ]; then
    echo "  clean.fa 已存在，跳过清洗: $CLEAN" | tee -a $LOG
else
    awk '/^>/{print $1; next} {print}' "$RAW" > "$CLEAN.tmp"
    rl=$(wc -l < "$RAW"); tl=$(wc -l < "$CLEAN.tmp")
    [ "$rl" != "$tl" ] && { echo "  行数不符 raw=$rl tmp=$tl，放弃" | tee -a $LOG; \rm -f "$CLEAN.tmp"; exit 4; }
    head -1 "$CLEAN.tmp" | grep -q ' ' && { echo "  header 仍含空格，放弃" | tee -a $LOG; \rm -f "$CLEAN.tmp"; exit 5; }
    mv "$CLEAN.tmp" "$CLEAN"
fi
echo "  clean.fa: $(du -h $CLEAN | cut -f1)  contigs=$(grep -c '^>' $CLEAN)  首行=$(head -1 $CLEAN)" | tee -a $LOG
echo "  软屏蔽比: $(awk '!/^>/{t+=length($0); n+=gsub(/[acgtn]/,"")} END{printf "%.1f%%", 100*n/t}' $CLEAN)" | tee -a $LOG

# ---------------- 4. 追加进物种清单 ----------------
if grep -qx "$SP" "$SPECIES_LIST_FILE"; then
    echo "  已在物种清单中" | tee -a $LOG
else
    cp "$SPECIES_LIST_FILE" "${SPECIES_LIST_FILE}.bak_$(date +%s)"
    echo "$SP" >> "$SPECIES_LIST_FILE"
    echo "  已追加进物种清单，现共 $(wc -l < $SPECIES_LIST_FILE) 株" | tee -a $LOG
fi

# ---------------- 5. 预检 ----------------
echo "" | tee -a $LOG
[ -s "$PROT_DB" ]           && echo "  蛋白库 OK: $PROT_DB ($(du -h $PROT_DB | cut -f1))" | tee -a $LOG || { echo "  蛋白库缺失"; exit 6; }
[ -s "$BRAKER_SCRIPT" ]     && echo "  braker.pl OK" | tee -a $LOG || { echo "  braker.pl 缺失"; exit 6; }
[ -d "$BUSCO_LINEAGE_DIR" ] && echo "  BUSCO 库 OK: $BUSCO_LINEAGE_DIR" | tee -a $LOG || { echo "  BUSCO 库缺失"; exit 6; }
echo "  可用空间: $(df -h "$DATA_ROOT" | awk 'NR==2{print $4}')  (需 >30G)" | tee -a $LOG

# ---------------- 6. 启动 ----------------
echo "" | tee -a $LOG
echo "[$(ts)] 启动 BRAKER2（8 线程）" | tee -a $LOG
nohup bash "$here/run_braker2_one.sh" "$SP" >> "$LOG_DIR/${SP}.launch.out" 2>&1 &
echo "  pid=$!" | tee -a $LOG
sleep 20
echo "  state: $(cat $STATE_DIR/${SP}.state 2>/dev/null || echo '尚未写入')" | tee -a $LOG
echo "" | tee -a $LOG
echo "  监控: watch -n 30 -t 'bash $here/monitor_addon_two.sh'   (先把脚本里 SPECIES 改成 $SP)" | tee -a $LOG
echo "  或:   tail -f $BATCH_ROOT/$SP/braker.log" | tee -a $LOG
