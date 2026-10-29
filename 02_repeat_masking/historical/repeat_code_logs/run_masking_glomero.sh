#!/bin/bash

# ================= 环境变量与路径 =================
export PATH="/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatModeler:/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatMasker:$PATH"

BASE_DIR="/media/desk16/tly8529/mucoromycota_annotation_braker"
INPUT_DIR="$BASE_DIR/01_Glomeromycotina_25"
OUT_MASKED="$BASE_DIR/03_softmasked_genomes"
OUT_RM_WORK="$BASE_DIR/04_repeat_working_dir"
LOG_DIR="$BASE_DIR/05_parallel_logs"

# 确保所有目录存在
mkdir -p "$OUT_MASKED" "$OUT_RM_WORK" "$LOG_DIR"

# ================= 单株处理函数 =================
process_genome() {
    GENOME="$1"
    FILENAME=$(basename "$GENOME")
    PREFIX="${FILENAME%.*}"
    LOG_FILE="$LOG_DIR/${PREFIX}.log"

    # 防重复与断点续跑：如果最终文件已存在，直接跳过
    if [ -f "$OUT_MASKED/${PREFIX}_softmasked.fa" ]; then
        echo "[跳过] $PREFIX 已存在软屏蔽结果。"
        return 0
    fi

    # 初始化日志
    echo "=====================================================" > "$LOG_FILE"
    echo "开始处理: $PREFIX" >> "$LOG_FILE"
    date >> "$LOG_FILE"

    # 进入专属工作目录
    WORK_DIR="$OUT_RM_WORK/$PREFIX"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR" || return 1

    echo "[1/3] 执行 BuildDatabase..." >> "$LOG_FILE"
    BuildDatabase -name "$PREFIX" "$GENOME" >> "$LOG_FILE" 2>&1

    echo "[2/3] 执行 RepeatModeler..." >> "$LOG_FILE"
    RepeatModeler -database "$PREFIX" -pa 6 >> "$LOG_FILE" 2>&1

    # 验证 families.fa 是否成功生成
    if [ -f "${PREFIX}-families.fa" ]; then
        echo "[3/3] 执行 RepeatMasker 软屏蔽 (-xsmall)..." >> "$LOG_FILE"
        RepeatMasker -lib "${PREFIX}-families.fa" -pa 6 -xsmall -gff -dir "$WORK_DIR" "$GENOME" >> "$LOG_FILE" 2>&1

        # 验证 .masked 文件是否成功生成并复制到总目录
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

# 导出函数和变量以供 xargs 并行调用
export -f process_genome
export OUT_MASKED OUT_RM_WORK LOG_DIR

# ================= 并行执行区 =================
echo "开始搜寻 $INPUT_DIR 下的 .gna 文件并提交并行处理..."

# 查找目标文件夹下的基因组（优先匹配 .gna，同时兼容 .fa/.fasta/.fna），送入并行队列
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.gna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | \
xargs -I {} -P 5 bash -c 'process_genome "$@"' _ {}

echo "Glomeromycotina 25株任务分配完成！"
