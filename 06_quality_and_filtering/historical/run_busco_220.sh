#!/usr/bin/env bash
set -u

source ~/bioinfo/miniforge3/bin/activate
conda activate busco6

GENOME_DIR=~/Desktop/220gna
PROJECT=~/bioinfo/mucoromycota_busco_local
DL_DIR=$PROJECT/downloads
LOG_DIR=$PROJECT/logs
OUT1=$PROJECT/results_mucoromycota
OUT2=$PROJECT/results_fungi
THREADS=16
MASTER_LOG=$LOG_DIR/run_busco_220.master.log

mkdir -p "$DL_DIR" "$LOG_DIR" "$OUT1" "$OUT2"

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

        if [ -d "$outdir/$sample" ]; then
            echo "[$(date)] SKIP $sample $lineage (result exists)" | tee -a "$MASTER_LOG"
            continue
        fi

        echo "[$(date)] START $sample $lineage" | tee -a "$MASTER_LOG"

        busco -i "$f" \
          -m genome \
          -l "$lineage" \
          -o "$sample" \
          --download_path "$DL_DIR" \
          --out_path "$outdir" \
          -c "$THREADS" \
          > "$LOG_DIR/${sample}.${lineage}.log" 2>&1

        status=$?
        if [ "$status" -ne 0 ]; then
            echo "[$(date)] FAIL  $sample $lineage (exit=$status)" | tee -a "$MASTER_LOG"
        else
            echo "[$(date)] END   $sample $lineage" | tee -a "$MASTER_LOG"
        fi
    done

    echo "===== END $lineage =====" | tee -a "$MASTER_LOG"
}

run_lineage "mucoromycota_odb12" "$OUT1"
run_lineage "fungi_odb12" "$OUT2"

echo "[$(date)] ALL BUSCO RUNS FINISHED" | tee -a "$MASTER_LOG"
