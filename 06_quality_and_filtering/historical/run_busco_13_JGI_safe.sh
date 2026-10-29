#!/usr/bin/env bash
set -u

source ~/bioinfo/miniforge3/bin/activate
conda activate busco6

GENOME_DIR=~/Desktop/13_JGI
RUN_TAG=$(date +%Y%m%d_%H%M%S)
PROJECT=~/bioinfo/mucoromycota_busco_local_13_JGI_$RUN_TAG
DL_DIR=~/Desktop/busco_db_mac
LOG_DIR=$PROJECT/logs
OUT1=$PROJECT/results_mucoromycota
OUT2=$PROJECT/results_fungi
THREADS=16
MASTER_LOG=$LOG_DIR/run_busco_13_JGI_offline.master.log

mkdir -p "$LOG_DIR" "$OUT1" "$OUT2"

cleanup() {
  echo "[$(date)] INTERRUPTED: stopping BUSCO batch" | tee -a "$MASTER_LOG"
  exit 130
}
trap cleanup INT TERM

run_lineage () {
  local lineage="$1"
  local outdir="$2"

  echo "===== START $lineage =====" | tee -a "$MASTER_LOG"

  for f in "$GENOME_DIR"/*.fna; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    sample=${base%.fna}
    sample=${sample%_genomic}

    if [ -d "$outdir/$sample" ] && ls "$outdir/$sample"/short_summary* >/dev/null 2>&1; then
      echo "[$(date)] SKIP $sample $lineage (complete result exists)" | tee -a "$MASTER_LOG"
      continue
    fi

    if [ -d "$outdir/$sample" ]; then
      mv "$outdir/$sample" "$outdir/${sample}.old_$(date +%Y%m%d_%H%M%S)"
    fi

    echo "[$(date)] START $sample $lineage" | tee -a "$MASTER_LOG"

    busco -i "$f" \
      -m genome \
      -l "$lineage" \
      --offline \
      --download_path "$DL_DIR" \
      -o "$sample" \
      --out_path "$outdir" \
      -c "$THREADS" \
      > "$LOG_DIR/${sample}.${lineage}.log" 2>&1

    status=$?
    if [ "$status" -ne 0 ]; then
      echo "[$(date)] FAIL $sample $lineage (exit=$status)" | tee -a "$MASTER_LOG"
    else
      echo "[$(date)] END $sample $lineage" | tee -a "$MASTER_LOG"
    fi
  done

  echo "===== END $lineage =====" | tee -a "$MASTER_LOG"
}

run_lineage "mucoromycota_odb12" "$OUT1"
run_lineage "fungi_odb12" "$OUT2"

echo "[$(date)] ALL BUSCO RUNS FINISHED" | tee -a "$MASTER_LOG"
