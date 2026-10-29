#!/bin/bash

# ================= 环境变量与路径 =================
export PATH="/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatModeler:/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatMasker:$PATH"

BASE_DIR="/media/desk16/tly8529/mucoromycota_annotation_braker"

# 要跑的两个输入目录（67株 + 119株 = 186株）
INPUT_DIR_1="$BASE_DIR/01_Mortierellomycotina_67"
INPUT_DIR_2="$BASE_DIR/01_Mucoromycotina_119"

# ================= 全新隔离的输出目录 =================
# 为了绝对保护原来的数据，这里建立全新的 v2 目录
OUT_MASKED_V2="$BASE_DIR/03_softmasked_genomes_v2_redo"
OUT_RM_WORK_V2="$BASE_DIR/04_repeat_working_dir_v2_redo"
LOG_DIR_V2="$BASE_DIR/05_parallel_logs_v2_redo"

mkdir -p "$OUT_MASKED_V2" "$OUT_RM_WORK_V2" "$LOG_DIR_V2"

# ================= 单株处理函数 =================
process_genome() {
    GENOME="$1"
    FILENAME=$(basename "$GENOME")
    PREFIX="${FILENAME%.*}"
    LOG_FILE="$LOG_DIR_V2/${PREFIX}.log"

    # 跳过在新目录中已经完成的基因组
    if [ -f "$OUT_MASKED_V2/${PREFIX}_softmasked.fa" ]; then
        echo "[跳过] $PREFIX 在重跑目录中已存在软屏蔽结果。"
        return 0
    fi

    echo "=====================================================" > "$LOG_FILE"
    echo "开始重跑处理: $PREFIX" >> "$LOG_FILE"
    date >> "$LOG_FILE"

    WORK_DIR="$OUT_RM_WORK_V2/$PREFIX"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR" || return 1

    echo "[1/3] 执行 BuildDatabase..." >> "$LOG_FILE"
    BuildDatabase -name "$PREFIX" "$GENOME" >> "$LOG_FILE" 2>&1

    # 核心参数：6线程
    echo "[2/3] 执行 RepeatModeler..." >> "$LOG_FILE"
    RepeatModeler -database "$PREFIX" -threads 6 >> "$LOG_FILE" 2>&1

    if [ -f "${PREFIX}-families.fa" ]; then
        # 核心参数：6线程，-xsmall
        echo "[3/3] 执行 RepeatMasker 软屏蔽 (-xsmall)..." >> "$LOG_FILE"
        RepeatMasker -lib "${PREFIX}-families.fa" -pa 6 -xsmall -gff -dir "$WORK_DIR" "$GENOME" >> "$LOG_FILE" 2>&1

        if [ -f "${FILENAME}.masked" ]; then
            # 将生成的软屏蔽文件放入 v2 隔离目录中
            cp "${FILENAME}.masked" "$OUT_MASKED_V2/${PREFIX}_softmasked.fa"
            echo "[成功] $PREFIX 重跑处理完成！" >> "$LOG_FILE"
        else
            echo "[错误] $PREFIX 未生成 .masked 文件，请检查日志。" >> "$LOG_FILE"
        fi
    else
        echo "[错误] $PREFIX 未生成 RepeatModeler 家族库，请检查日志。" >> "$LOG_FILE"
    fi
    date >> "$LOG_FILE"
}

export -f process_genome
export OUT_MASKED_V2 OUT_RM_WORK_V2 LOG_DIR_V2

# ================= 并行执行区 =================
echo "开始搜寻指定的两批基因组 (共 186 株) 并提交并行处理..."
# 核心参数：10个并行 (-P 10)
find "$INPUT_DIR_1" "$INPUT_DIR_2" -maxdepth 1 -type f \( -name "*.gna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | \
xargs -I {} -P 10 bash -c 'process_genome "$@"' _ {}

echo "186株基因组重跑任务分配完成！"
