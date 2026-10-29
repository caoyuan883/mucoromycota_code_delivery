#!/usr/bin/env bash

set -u
set -o pipefail

PROTEIN_DIR="/home/py/Desktop/Mucoromycota_BRAKER_longest_isoform/02_longest_isoform_proteins"

PROJECT="$HOME/bioinfo/mucoromycota_busco_longest_isoform_test3"

DL_DIR="$HOME/Desktop/busco_db_mac"

LOG_DIR="$PROJECT/logs"

SAMPLES=(
    "Cokeromyces_recurvatus"
    "Cryptocalciella_humilis"
    "Cunninghamella_bertholletiae"
)

LINEAGES=(
    "mucoromycota_odb12"
    "fungi_odb12"
)

mkdir -p \
    "$PROJECT/results_mucoromycota" \
    "$PROJECT/results_fungi" \
    "$LOG_DIR"

echo "============================================================"
echo "BUSCO longest-isoform测试"
echo "时间：$(date)"
echo "输入目录：$PROTEIN_DIR"
echo "输出目录：$PROJECT"
echo "BUSCO数据库：$DL_DIR"
echo "============================================================"

if ! command -v busco >/dev/null 2>&1; then
    echo "错误：当前环境找不到busco命令"
    echo "请先激活此前使用的BUSCO 6环境"
    exit 1
fi

echo "BUSCO版本：$(busco --version 2>&1 | head -1)"

failed=0
completed=0

for sample in "${SAMPLES[@]}"; do

    input="$PROTEIN_DIR/${sample}_longest_isoform.faa"

    if [ ! -s "$input" ]; then
        echo "ERROR：输入文件不存在或为空：$input"
        failed=$((failed + 2))
        continue
    fi

    echo
    echo "样本：$sample"
    echo "蛋白数：$(grep -c '^>' "$input")"

    for lineage in "${LINEAGES[@]}"; do

        if [ "$lineage" = "mucoromycota_odb12" ]; then
            out_path="$PROJECT/results_mucoromycota"
        else
            out_path="$PROJECT/results_fungi"
        fi

        out_name="${sample}_longest_isoform"

        existing_summary=$(find \
            "$out_path/$out_name" \
            -maxdepth 1 \
            -type f \
            -name 'short_summary*.txt' \
            -print -quit 2>/dev/null || true)

        if [ -n "$existing_summary" ]; then
            echo "SKIP：已有结果 $sample $lineage"
            completed=$((completed + 1))
            continue
        fi

        if [ -d "$out_path/$out_name" ]; then
            rm -rf "$out_path/$out_name"
        fi

        log="$LOG_DIR/${sample}.${lineage}.log"

        echo "[$(date)] START $sample $lineage"

        busco \
            -i "$input" \
            -m proteins \
            -l "$lineage" \
            --offline \
            --download_path "$DL_DIR" \
            -o "$out_name" \
            --out_path "$out_path" \
            -c 16 \
            > "$log" 2>&1

        status=$?

        if [ "$status" -eq 0 ]; then
            echo "[$(date)] END   $sample $lineage"
            completed=$((completed + 1))
        else
            echo "[$(date)] FAIL  $sample $lineage exit=$status"
            echo "日志：$log"
            failed=$((failed + 1))
        fi

    done
done

echo
echo "============================================================"
echo "全部任务结束：$(date)"
echo "成功/已有：$completed"
echo "失败：$failed"
echo "预期总任务：6"
echo "============================================================"

if [ "$failed" -ne 0 ]; then
    exit 1
fi
