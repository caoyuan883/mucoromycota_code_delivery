#!/usr/bin/env bash

set -u
set -o pipefail

PROTEIN_DIR="/home/py/Desktop/Mucoromycota_BRAKER_longest_isoform/02_longest_isoform_proteins"

PROJECT="$HOME/bioinfo/mucoromycota_busco_longest_isoform_213"

DL_DIR="$HOME/Desktop/busco_db_mac"

LOG_DIR="$PROJECT/logs"

STATUS_FILE="$PROJECT/busco_213_status.tsv"

THREADS=16

LINEAGES=(
    "mucoromycota_odb12"
    "fungi_odb12"
)

mkdir -p \
    "$PROJECT/results_mucoromycota" \
    "$PROJECT/results_fungi" \
    "$LOG_DIR"

echo -e \
"sample\tlineage\tstatus\texit_code\tinput\toutput\tlog" \
> "$STATUS_FILE"

echo "================================================================"
echo "BUSCO longest-isoform 213株批处理"
echo "开始时间：$(date)"
echo "输入目录：$PROTEIN_DIR"
echo "输出目录：$PROJECT"
echo "数据库目录：$DL_DIR"
echo "模式：proteins"
echo "核心数：$THREADS"
echo "数据库：${LINEAGES[*]}"
echo "================================================================"

if ! command -v busco >/dev/null 2>&1; then
    echo "错误：当前环境找不到busco命令"
    echo "请先激活此前使用的BUSCO 6环境"
    exit 1
fi

echo "BUSCO版本：$(busco --version 2>&1 | head -1)"

if [ ! -d "$PROTEIN_DIR" ]; then
    echo "错误：输入目录不存在：$PROTEIN_DIR"
    exit 1
fi

if [ ! -d "$DL_DIR" ]; then
    echo "错误：BUSCO数据库目录不存在：$DL_DIR"
    exit 1
fi

for lineage in "${LINEAGES[@]}"; do
    if [ ! -d "$DL_DIR/lineages/$lineage" ] &&
       [ ! -d "$DL_DIR/$lineage" ]; then
        echo "错误：找不到BUSCO数据库：$lineage"
        echo "已检查："
        echo "  $DL_DIR/lineages/$lineage"
        echo "  $DL_DIR/$lineage"
        exit 1
    fi
done

mapfile -d '' PROTEIN_FILES < <(
    find "$PROTEIN_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*_longest_isoform.faa' \
        -print0 |
    sort -z
)

FILE_COUNT="${#PROTEIN_FILES[@]}"

echo "检测到longest-isoform蛋白文件：$FILE_COUNT"

if [ "$FILE_COUNT" -ne 213 ]; then
    echo "错误：预期213个蛋白文件，实际检测到$FILE_COUNT"
    echo "为防止不完整批处理，脚本停止。"
    exit 1
fi

empty_count=0

for input in "${PROTEIN_FILES[@]}"; do
    if [ ! -s "$input" ]; then
        echo "空文件：$input"
        empty_count=$((empty_count + 1))
    fi
done

if [ "$empty_count" -ne 0 ]; then
    echo "错误：发现$empty_count个空蛋白文件"
    exit 1
fi

total_tasks=$((FILE_COUNT * ${#LINEAGES[@]}))
completed=0
skipped=0
failed=0
task_number=0

for input in "${PROTEIN_FILES[@]}"; do

    filename=$(basename "$input")
    sample="${filename%_longest_isoform.faa}"

    protein_count=$(grep -c '^>' "$input" || true)

    echo
    echo "----------------------------------------------------------------"
    echo "样本：$sample"
    echo "蛋白数：$protein_count"
    echo "----------------------------------------------------------------"

    if [ "$protein_count" -eq 0 ]; then
        echo "ERROR：没有检测到FASTA记录：$input"

        for lineage in "${LINEAGES[@]}"; do
            echo -e \
"${sample}\t${lineage}\tFAILED_NO_FASTA_RECORDS\t1\t${input}\t\t" \
            >> "$STATUS_FILE"

            failed=$((failed + 1))
        done

        continue
    fi

    for lineage in "${LINEAGES[@]}"; do

        task_number=$((task_number + 1))

        if [ "$lineage" = "mucoromycota_odb12" ]; then
            out_path="$PROJECT/results_mucoromycota"
        else
            out_path="$PROJECT/results_fungi"
        fi

        out_name="${sample}_longest_isoform"
        result_dir="$out_path/$out_name"
        log="$LOG_DIR/${sample}.${lineage}.log"

        existing_summary=$(find \
            "$result_dir" \
            -maxdepth 1 \
            -type f \
            -name 'short_summary*.txt' \
            -print -quit 2>/dev/null || true)

        if [ -n "$existing_summary" ]; then
            if grep -qE \
                'C:[[:space:]]*[0-9.]+%' \
                "$existing_summary"; then

                skipped=$((skipped + 1))

                echo \
"[$task_number/$total_tasks] SKIP $sample $lineage"

                echo -e \
"${sample}\t${lineage}\tSKIPPED_COMPLETE\t0\t${input}\t${result_dir}\t${log}" \
                >> "$STATUS_FILE"

                continue
            fi
        fi

        if [ -d "$result_dir" ]; then
            echo \
"删除不完整旧结果：$result_dir"

            rm -rf "$result_dir"
        fi

        echo \
"[$task_number/$total_tasks] START $sample $lineage $(date)"

        busco \
            -i "$input" \
            -m proteins \
            -l "$lineage" \
            --offline \
            --download_path "$DL_DIR" \
            -o "$out_name" \
            --out_path "$out_path" \
            -c "$THREADS" \
            > "$log" 2>&1

        status=$?

        if [ "$status" -eq 0 ]; then

            summary=$(find \
                "$result_dir" \
                -maxdepth 1 \
                -type f \
                -name 'short_summary*.txt' \
                -print -quit 2>/dev/null || true)

            if [ -n "$summary" ] &&
               grep -qE \
                   'C:[[:space:]]*[0-9.]+%' \
                   "$summary"; then

                completed=$((completed + 1))

                echo \
"[$task_number/$total_tasks] END $sample $lineage $(date)"

                echo -e \
"${sample}\t${lineage}\tCOMPLETED\t0\t${input}\t${result_dir}\t${log}" \
                >> "$STATUS_FILE"

            else
                failed=$((failed + 1))

                echo \
"[$task_number/$total_tasks] FAIL_NO_SUMMARY $sample $lineage"

                echo -e \
"${sample}\t${lineage}\tFAILED_NO_VALID_SUMMARY\t2\t${input}\t${result_dir}\t${log}" \
                >> "$STATUS_FILE"
            fi

        else
            failed=$((failed + 1))

            echo \
"[$task_number/$total_tasks] FAIL $sample $lineage exit=$status"

            echo "日志：$log"

            echo -e \
"${sample}\t${lineage}\tFAILED_BUSCO\t${status}\t${input}\t${result_dir}\t${log}" \
            >> "$STATUS_FILE"
        fi

    done
done

valid_mu=$(find \
    "$PROJECT/results_mucoromycota" \
    -type f \
    -name 'short_summary*.txt' |
wc -l)

valid_fu=$(find \
    "$PROJECT/results_fungi" \
    -type f \
    -name 'short_summary*.txt' |
wc -l)

echo
echo "================================================================"
echo "BUSCO批处理结束"
echo "结束时间：$(date)"
echo "输入蛋白文件：$FILE_COUNT"
echo "总任务数：$total_tasks"
echo "本次完成：$completed"
echo "已有跳过：$skipped"
echo "失败：$failed"
echo "Mucoromycota summary：$valid_mu"
echo "Fungi summary：$valid_fu"
echo "状态表：$STATUS_FILE"
echo "================================================================"

if [ "$failed" -ne 0 ]; then
    exit 1
fi
