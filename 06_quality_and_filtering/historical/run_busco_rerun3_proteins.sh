source ~/bioinfo/miniforge3/bin/activate && conda activate busco6

PROTEIN_DIR="/home/py/Desktop/3_rerun_protein"
PROJECT="$HOME/bioinfo/mucoromycota_busco_protein_rerun3"
DL_DIR="$HOME/Desktop/busco_db_mac"
LOG_DIR="$PROJECT/logs"

mkdir -p \
"$LOG_DIR" \
"$PROJECT/results_mucoromycota" \
"$PROJECT/results_fungi"

echo "发现蛋白文件：$(find "$PROTEIN_DIR" -maxdepth 1 -type f -name '*_protein.faa' | wc -l)"

for f in "$PROTEIN_DIR"/*_protein.faa; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    sample="${base%_protein.faa}"

    for lineage in mucoromycota_odb12 fungi_odb12; do

        if [ "$lineage" = "mucoromycota_odb12" ]; then
            outdir="$PROJECT/results_mucoromycota"
        else
            outdir="$PROJECT/results_fungi"
        fi

        if [ -d "$outdir/$sample" ] &&
           find "$outdir/$sample" -maxdepth 1 \
             -name 'short_summary*.txt' -print -quit |
             grep -q .; then
            echo "[$(date)] SKIP $sample $lineage"
            continue
        fi

        [ -d "$outdir/$sample" ] &&
            rm -rf "$outdir/$sample"

        echo "[$(date)] START $sample $lineage"

        busco \
          -i "$f" \
          -m proteins \
          -l "$lineage" \
          --offline \
          --download_path "$DL_DIR" \
          -o "$sample" \
          --out_path "$outdir" \
          -c 16 \
          > "$LOG_DIR/${sample}.${lineage}.log" 2>&1

        status=$?

        if [ "$status" -eq 0 ]; then
            echo "[$(date)] END $sample $lineage"
        else
            echo "[$(date)] FAIL $sample $lineage exit=$status"
            echo "日志：$LOG_DIR/${sample}.${lineage}.log"
        fi
    done
done

echo "[$(date)] ALL BUSCO RUNS FINISHED"
