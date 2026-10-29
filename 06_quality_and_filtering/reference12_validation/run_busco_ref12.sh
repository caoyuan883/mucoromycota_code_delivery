#!/usr/bin/env bash
# [release refactor] activate the BUSCO environment yourself before running.
# historical: source ~/bioinfo/miniforge3/bin/activate && conda activate busco6
# historical: IN_DIR=/home/py/Desktop/12_faa_12_gff_technical_validation/12_faa_longest_isoform
# historical: PROJECT=$HOME/bioinfo/benchmark12_busco_ref_protein ; DL_DIR=$HOME/Desktop/busco_db_mac
IN_DIR="${REF12_BASE:?set REF12_BASE}/12_faa_longest_isoform"
PROJECT="${REF12_BUSCO_PROJECT:?set REF12_BUSCO_PROJECT}"
DL_DIR="${BUSCO_DOWNLOAD_DIR:?set BUSCO_DOWNLOAD_DIR}"
LOG_DIR="$PROJECT/logs"
OUT="$PROJECT/results_mucoromycota"
mkdir -p "$LOG_DIR" "$OUT"

echo "发现蛋白文件：$(ls "$IN_DIR"/*_longest_isoform.faa 2>/dev/null | wc -l)"

for f in "$IN_DIR"/*_longest_isoform.faa; do
  [ -e "$f" ] || continue
  base=$(basename "$f"); sample="${base%_longest_isoform.faa}"

  if [ -d "$OUT/$sample" ] && \
     find "$OUT/$sample" -maxdepth 1 -name 'short_summary*.txt' -print -quit | grep -q .; then
    echo "[$(date)] SKIP $sample"; continue
  fi
  [ -d "$OUT/$sample" ] && rm -rf "$OUT/$sample"

  echo "[$(date)] START $sample"
  busco -i "$f" -m proteins -l mucoromycota_odb12 \
        --offline --download_path "$DL_DIR" \
        -o "$sample" --out_path "$OUT" -c 16 \
        > "$LOG_DIR/${sample}.log" 2>&1
  [ $? -eq 0 ] && echo "[$(date)] END $sample" \
                || echo "[$(date)] FAIL $sample  日志：$LOG_DIR/${sample}.log"
done
echo "[$(date)] ALL DONE"
