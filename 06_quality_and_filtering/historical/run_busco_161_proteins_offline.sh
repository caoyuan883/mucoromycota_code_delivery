#!/usr/bin/env bash
set -u

source ~/bioinfo/miniforge3/bin/activate
conda activate busco6

PROT_DIR=~/Desktop/161_proteins
PROJECT=~/bioinfo/protein_busco_161
DL_DIR=~/Desktop/busco_db_mac
LOG_DIR=$PROJECT/logs
OUT1=$PROJECT/results_mucoromycota
OUT2=$PROJECT/results_fungi
THREADS=16
MASTER_LOG=$LOG_DIR/run_busco_161_proteins_offline.master.log

mkdir -p "$LOG_DIR" "$OUT1" "$OUT2"

cleanup() {
    echo "[$(date)] INTERRUPTED: stopping BUSCO protein batch" | tee -a "$MASTER_LOG"
    exit 130
}
trap cleanup INT TERM

run_lineage () {
    local lineage="$1"
    local outdir="$2"

    echo "===== START $lineage =====" | tee -a "$MASTER_LOG"

    for f in "$PROT_DIR"/*.faa; do
        [ -e "$f" ] || continue
        base=$(basename "$f")
        sample=${base%.faa}
        sample=${sample%_protein}

        if [ -d "$outdir/$sample" ] && ls "$outdir/$sample"/short_summary* >/dev/null 2>&1; then
            echo "[$(date)] SKIP $sample $lineage (complete result exists)" | tee -a "$MASTER_LOG"
            continue
        fi

        [ -d "$outdir/$sample" ] && rm -rf "$outdir/$sample"

        echo "[$(date)] START $sample $lineage" | tee -a "$MASTER_LOG"

        busco -i "$f" \
          -m proteins \
          -l "$lineage" \
          --offline \
          --download_path "$DL_DIR" \
          -o "$sample" \
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

echo "[$(date)] ALL PROTEIN BUSCO RUNS FINISHED" | tee -a "$MASTER_LOG"
