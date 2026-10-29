#!/usr/bin/env bash
source ~/bioinfo/miniforge3/bin/activate && conda activate busco6

BASE="/home/py/Desktop/5_outgroup"
FAA_DIR="$BASE/5_faa"
OUT="$BASE/busco_protein"
LOG_DIR="$BASE/logs"
DL_DIR="$HOME/Desktop/busco_db_mac"
LINEAGE="fungi_odb12"

mkdir -p "$OUT" "$LOG_DIR"
echo "BUSCO: $(busco --version 2>&1 | head -1)"
echo "发现蛋白文件：$(ls "$FAA_DIR"/*_protein.faa 2>/dev/null | wc -l)"
echo

for f in "$FAA_DIR"/*_protein.faa; do
  [ -e "$f" ] || continue
  base=$(basename "$f"); sample="${base%_protein.faa}"

  if [ -d "$OUT/$sample" ] && \
     find "$OUT/$sample" -maxdepth 1 -name 'short_summary*.txt' -print -quit | grep -q .; then
    echo "[$(date +%H:%M:%S)] SKIP $sample"; continue
  fi
  [ -d "$OUT/$sample" ] && rm -rf "$OUT/$sample"

  echo "[$(date +%H:%M:%S)] START $sample  ($(grep -c '^>' "$f") 条蛋白)"
  busco \
    -i "$f" \
    -m proteins \
    -l "$LINEAGE" \
    --offline \
    --download_path "$DL_DIR" \
    -o "$sample" \
    --out_path "$OUT" \
    -c 16 \
    > "$LOG_DIR/${sample}.protein.${LINEAGE}.log" 2>&1

  [ $? -eq 0 ] && echo "[$(date +%H:%M:%S)] END   $sample" \
               || echo "[$(date +%H:%M:%S)] FAIL  $sample  见 $LOG_DIR/${sample}.protein.${LINEAGE}.log"
done
echo; echo "[$(date +%H:%M:%S)] 全部结束"
