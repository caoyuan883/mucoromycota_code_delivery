#!/usr/bin/env bash
set -u

source ~/bioinfo/miniforge3/bin/activate
conda activate busco6

PROTEIN_DIR="/home/py/Desktop/Mucoromycota_BRAKER_core_proteomes"
PROJECT="$HOME/bioinfo/mucoromycota_busco_protein_local"
DL_DIR="$HOME/Desktop/busco_db_mac"

LOG_DIR="$PROJECT/logs"
OUT1="$PROJECT/results_mucoromycota"
OUT2="$PROJECT/results_fungi"

THREADS=16
MASTER_LOG="$LOG_DIR/run_busco_proteins_offline.master.log"

mkdir -p "$LOG_DIR" "$OUT1" "$OUT2"

if [ ! -d "$PROTEIN_DIR" ]; then
    echo "ERROR: 找不到蛋白质目录：$PROTEIN_DIR"
    exit 1
fi

if [ ! -d "$DL_DIR" ]; then
    echo "ERROR: 找不到 BUSCO 离线数据库：$DL_DIR"
    exit 1
fi

N=$(find "$PROTEIN_DIR" -maxdepth 1 -type f \
    -name '*_protein.faa' | wc -l)

echo "发现蛋白质文件：$N" | tee -a "$MASTER_LOG"

cleanup() {
    echo "[$(date)] INTERRUPTED: stopping BUSCO batch" |
        tee -a "$MASTER_LOG"
    exit 130
}

trap cleanup INT TERM

run_lineage() {
    lineage="$1"
    outdir="$2"

    echo "===== START $lineage =====" | tee -a "$MASTER_LOG"

    while IFS= read -r -d '' f; do
        base=$(basename "$f")
        sample="${base%_protein.faa}"

        sample_dir="$outdir/$sample"
        sample_log="$LOG_DIR/${sample}.${lineage}.log"

        if [ -d "$sample_dir" ] &&
           find "$sample_dir" -maxdepth 1 \
             -name 'short_summary*.txt' -print -quit |
             grep -q .; then

            echo "[$(date)] SKIP $sample $lineage" |
                tee -a "$MASTER_LOG"
            continue
        fi

        if [ -d "$sample_dir" ]; then
            echo "[$(date)] REMOVE incomplete: $sample $lineage" |
                tee -a "$MASTER_LOG"
            rm -rf "$sample_dir"
        fi

        echo "[$(date)] START $sample $lineage" |
            tee -a "$MASTER_LOG"

        busco \
            -i "$f" \
            -m proteins \
            -l "$lineage" \
            --offline \
            --download_path "$DL_DIR" \
            -o "$sample" \
            --out_path "$outdir" \
            -c "$THREADS" \
            > "$sample_log" 2>&1

        status=$?

        if [ "$status" -eq 0 ]; then
            echo "[$(date)] END $sample $lineage" |
                tee -a "$MASTER_LOG"
        else
            echo "[$(date)] FAIL $sample $lineage exit=$status" |
                tee -a "$MASTER_LOG"
        fi

    done < <(
        find "$PROTEIN_DIR" -maxdepth 1 -type f \
            -name '*_protein.faa' -print0 | sort -z
    )

    echo "===== END $lineage =====" | tee -a "$MASTER_LOG"
}

run_lineage "mucoromycota_odb12" "$OUT1"
run_lineage "fungi_odb12" "$OUT2"

echo "[$(date)] ALL BUSCO PROTEIN RUNS FINISHED" |
    tee -a "$MASTER_LOG"

echo "Mucoromycota 完成结果数：$(find "$OUT1" -name 'short_summary*.txt' | wc -l)"
echo "Fungi 完成结果数：$(find "$OUT2" -name 'short_summary*.txt' | wc -l)"
