#!/usr/bin/env bash
source ~/bioinfo/miniforge3/bin/activate && conda activate busco6

BASE="/home/py/Desktop/5_outgroup"
GENOME_DIR="$BASE/5_genome"
OUT="$BASE/busco_genome"
LOG_DIR="$BASE/logs"
DL_DIR="$HOME/Desktop/busco_db_mac"
LINEAGE="fungi_odb12"

mkdir -p "$OUT" "$LOG_DIR"
echo "BUSCO: $(busco --version 2>&1 | head -1)"
echo "发现基因组文件：$(ls "$GENOME_DIR"/*_genomic.fna 2>/dev/null | wc -l)"
echo

for f in "$GENOME_DIR"/*_genomic.fna; do
  [ -e "$f" ] || continue
  base=$(basename "$f"); sample="${base%_genomic.fna}"

  if [ -d "$OUT/$sample" ] && \
     find "$OUT/$sample" -maxdepth 1 -name 'short_summary*.txt' -print -quit | grep -q .; then
    echo "[$(date +%H:%M:%S)] SKIP $sample"; continue
  fi
  [ -d "$OUT/$sample" ] && rm -rf "$OUT/$sample"

  echo "[$(date +%H:%M:%S)] START $sample  ($(du -h "$f" | cut -f1))"
  busco \
    -i "$f" \
    -m genome \
    -l "$LINEAGE" \
    --offline \
    --download_path "$DL_DIR" \
    -o "$sample" \
    --out_path "$OUT" \
    -c 16 \
    > "$LOG_DIR/${sample}.genome.${LINEAGE}.log" 2>&1

  if [ $? -eq 0 ]; then
    echo "[$(date +%H:%M:%S)] END   $sample"
  else
    echo "[$(date +%H:%M:%S)] FAIL  $sample  见 $LOG_DIR/${sample}.genome.${LINEAGE}.log"
  fi
done

echo
echo "[$(date +%H:%M:%S)] 全部结束"
