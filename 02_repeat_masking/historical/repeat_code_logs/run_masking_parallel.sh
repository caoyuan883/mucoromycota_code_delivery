#!/bin/bash

export PATH="/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatModeler:/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatMasker:$PATH"

BASE_DIR="/media/desk16/tly8529/mucoromycota_annotation_braker"
OUT_MASKED="$BASE_DIR/03_softmasked_genomes"
OUT_RM_WORK="$BASE_DIR/04_repeat_working_dir"
LOG_DIR="$BASE_DIR/05_parallel_logs"
JOB_LIST="$BASE_DIR/genome_jobs.txt"

RMODELER_THREADS=6
RMASKER_PA=2
PARALLEL_JOBS=10

mkdir -p "$OUT_MASKED" "$OUT_RM_WORK" "$LOG_DIR"

: > "$JOB_LIST"

find "$BASE_DIR/01_Mucoromycotina_119" "$BASE_DIR/01_Mortierellomycotina_67" \
  -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | sort > "$JOB_LIST"

process_genome() {
    GENOME="$1"
    FILENAME=$(basename "$GENOME")
    PREFIX="${FILENAME%.*}"
    WORK_DIR="$OUT_RM_WORK/$PREFIX"
    LOG_FILE="$LOG_DIR/${PREFIX}.log"

    if [ -f "$OUT_MASKED/${PREFIX}_softmasked.fa" ]; then
        echo "[跳过] $PREFIX 已完成" >> "$LOG_FILE"
        return 0
    fi

    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR" || exit 1

    {
        echo "====================================================="
        echo "开始处理: $PREFIX"
        date
        echo "[1/3] BuildDatabase"
    } >> "$LOG_FILE"

    BuildDatabase -name "$PREFIX" "$GENOME" >> "$LOG_FILE" 2>&1

    {
        echo "[2/3] RepeatModeler"
        echo "RepeatModeler threads = $RMODELER_THREADS"
    } >> "$LOG_FILE"

    RepeatModeler -database "$PREFIX" -threads $RMODELER_THREADS >> "$LOG_FILE" 2>&1

    if [ -f "${PREFIX}-families.fa" ]; then
        {
            echo "[3/3] RepeatMasker"
            echo "RepeatMasker -pa = $RMASKER_PA"
        } >> "$LOG_FILE"

        RepeatMasker -lib "${PREFIX}-families.fa" -pa $RMASKER_PA -xsmall -gff -dir "$WORK_DIR" "$GENOME" >> "$LOG_FILE" 2>&1

        if [ -f "$WORK_DIR/${FILENAME}.masked" ]; then
            cp "$WORK_DIR/${FILENAME}.masked" "$OUT_MASKED/${PREFIX}_softmasked.fa"
            echo "[成功] $PREFIX 完成" >> "$LOG_FILE"
        else
            echo "[错误] $PREFIX 没生成 masked 文件" >> "$LOG_FILE"
        fi
    else
        echo "[错误] $PREFIX 没生成 families.fa" >> "$LOG_FILE"
    fi
}

export -f process_genome
export BASE_DIR OUT_MASKED OUT_RM_WORK LOG_DIR RMODELER_THREADS RMASKER_PA

cat "$JOB_LIST" | xargs -I {} -P $PARALLEL_JOBS bash -c 'process_genome "$@"' _ {}
