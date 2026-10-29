#!/bin/bash

# ================= 环境变量与路径 =================
export PATH="/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatModeler:/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatMasker:$PATH"

BASE_DIR="/media/desk16/tly8529/mucoromycota_annotation_braker"

# 最后的 7 个基因组输入目录
INPUT_DIR="$BASE_DIR/7_genomic"

# ================= 专属的输出目录 =================
OUT_MASKED="$BASE_DIR/03_softmasked_genomes_7"
OUT_RM_WORK="$BASE_DIR/04_repeat_working_dir_7"
LOG_DIR="$BASE_DIR/05_parallel_logs_7"

mkdir -p "$OUT_MASKED" "$OUT_RM_WORK" "$LOG_DIR"

# ================= 单株处理函数 =================
process_genome() {
    GENOME="$1"
    FILENAME=$(basename "$GENOME")
    PREFIX="${FILENAME%.*}"
    LOG_FILE="$LOG_DIR/${PREFIX}.log"

    if [ -f "$OUT_MASKED/${PREFIX}_softmasked.fa" ]; then
        echo "[跳过] $PREFIX 已存在软屏蔽结果。"
        return 0
    fi

    echo "=====================================================" > "$LOG_FILE"
    echo "开始处理最后的7株: $PREFIX" >> "$LOG_FILE"
    date >> "$LOG_FILE"

    WORK_DIR="$OUT_RM_WORK/$PREFIX"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR" || return 1

    echo "[1/3] 执行 BuildDatabase..." >> "$LOG_FILE"
    BuildDatabase -name "$PREFIX" "$GENOME" >> "$LOG_FILE" 2>&1

    echo "[2/3] 执行 RepeatModeler..." >> "$LOG_FILE"
    RepeatModeler -database "$PREFIX" -threads 6 >> "$LOG_FILE" 2>&1

    if [ -f "${PREFIX}-families.fa" ]; then
        echo "[3/3] 执行 RepeatMasker 软屏蔽 (-xsmall)..." >> "$LOG_FILE"
        RepeatMasker -lib "${PREFIX}-families.fa" -pa 6 -xsmall -gff -dir "$WORK_DIR" "$GENOME" >> "$LOG_FILE" 2>&1

        if [ -f "${FILENAME}.masked" ]; then
            cp "${FILENAME}.masked" "$OUT_MASKED/${PREFIX}_softmasked.fa"
            echo "[成功] $PREFIX 处理完成！" >> "$LOG_FILE"
        else
            echo "[错误] $PREFIX 未生成 .masked 文件，请检查日志。" >> "$LOG_FILE"
        fi
    else
        echo "[错误] $PREFIX 未生成 RepeatModeler 家族库，请检查日志。" >> "$LOG_FILE"
    fi
    date >> "$LOG_FILE"
}

export -f process_genome
export OUT_MASKED OUT_RM_WORK LOG_DIR

# ================= 并行执行区 =================
echo "开始搜寻 7_genomic 目录下的基因组并提交并行处理..."
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.gna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | \
xargs -I {} -P 10 bash -c 'process_genome "$@"' _ {}

echo "7株基因组任务分配完成！"
