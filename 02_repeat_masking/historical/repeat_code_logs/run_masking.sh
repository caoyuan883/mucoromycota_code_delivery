#!/bin/bash

export PATH="/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatModeler:/media/desk16/tly8529/glo_braker_test/08_env/miniconda3/envs/glo_repeat/share/RepeatMasker:$PATH"

BASE_DIR="/media/desk16/tly8529/mucoromycota_annotation_braker"

INPUT_DIRS=(
    "$BASE_DIR/01_Mucoromycotina_119"
    "$BASE_DIR/01_Mortierellomycotina_67"
)

OUT_MASKED="$BASE_DIR/03_softmasked_genomes"
OUT_RM_WORK="$BASE_DIR/04_repeat_working_dir"

RMODELER_THREADS=48
RMASKER_PA=12

mkdir -p "$OUT_MASKED"
mkdir -p "$OUT_RM_WORK"

for DIR in "${INPUT_DIRS[@]}"; do
    echo "====================================================="
    echo "开始扫描目录: $DIR"

    shopt -s nullglob
    for GENOME in "$DIR"/*.{fa,fasta,fna}; do
        FILENAME=$(basename "$GENOME")
        PREFIX="${FILENAME%.*}"

        if [ -f "$OUT_MASKED/${PREFIX}_softmasked.fa" ]; then
            echo "[跳过] $PREFIX 软屏蔽已完成。"
            continue
        fi

        echo "-----------------------------------------------------"
        echo "正在处理基因组: $PREFIX"

        WORK_DIR="$OUT_RM_WORK/$PREFIX"
        mkdir -p "$WORK_DIR"
        cd "$WORK_DIR" || exit 1

        echo "[1/3] 执行 BuildDatabase..."
        BuildDatabase -name "$PREFIX" "$GENOME"

        echo "[2/3] 执行 RepeatModeler..."
        RepeatModeler -database "$PREFIX" -threads $RMODELER_THREADS

        if [ -f "${PREFIX}-families.fa" ]; then
            echo "[3/3] 执行 RepeatMasker 软屏蔽 (-xsmall)..."
            RepeatMasker -lib "${PREFIX}-families.fa" -pa $RMASKER_PA -xsmall -gff -dir "$WORK_DIR" "$GENOME"

            if [ -f "$WORK_DIR/${FILENAME}.masked" ]; then
                cp "$WORK_DIR/${FILENAME}.masked" "$OUT_MASKED/${PREFIX}_softmasked.fa"
                echo "[成功] $PREFIX 处理完成！"
            else
                echo "[错误] $PREFIX 未生成 .masked 文件，请检查日志。"
            fi
        else
            echo "[错误] $PREFIX 未生成 RepeatModeler 家族库，请检查日志。"
        fi
    done
    shopt -u nullglob
done

echo "====================================================="
echo "所有指定目录下的基因组预处理已完成！"
