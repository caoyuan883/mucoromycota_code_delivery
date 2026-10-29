#!/bin/bash
# ============================================================================
# Step 01 — 统一清洗 125 株基因组 FASTA header
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

MAPPING=$DATA_ROOT/species_mode_mapping.csv

mkdir -p $GENOME_DIR
mkdir -p $BATCH_ROOT
LOG=$BATCH_ROOT/_clean_headers.log
: > $LOG

echo "===== Step 01: 清洗 125 株基因组 header =====" | tee -a $LOG
echo "  输入: $RAW_GENOME_DIR" | tee -a $LOG
echo "  输出: $GENOME_DIR" | tee -a $LOG
echo "  mapping: $MAPPING" | tee -a $LOG
echo "  开始时间: $(ts)" | tee -a $LOG
echo "" | tee -a $LOG

awk -F',' '$NF=="BRAKER2" {print $1}' $MAPPING > $SPECIES_LIST_FILE
n_species=$(wc -l < $SPECIES_LIST_FILE)
echo "  从 mapping 提取 BRAKER2 物种: $n_species 株" | tee -a $LOG
echo "" | tee -a $LOG

if [ "$n_species" -ne 125 ]; then
    echo "  ⚠️  警告: 期望 125 株，实际 $n_species 株" | tee -a $LOG
    exit 2
fi

ok=0
skip=0
fail=0
missing=0
i=0

for sp in $(cat $SPECIES_LIST_FILE); do
    i=$((i + 1))
    raw=$RAW_GENOME_DIR/${sp}_genomic_softmasked.fa
    clean=$GENOME_DIR/${sp}_genomic_softmasked.clean.fa
    
    if [ ! -s "$raw" ]; then
        echo "  [$i/$n_species] ❌ $sp: 原始文件缺失 $raw" | tee -a $LOG
        missing=$((missing + 1))
        continue
    fi
    
    if [ -s "$clean" ]; then
        raw_lines=$(wc -l < $raw)
        clean_lines=$(wc -l < $clean)
        if [ "$raw_lines" = "$clean_lines" ]; then
            echo "  [$i/$n_species] ⏭  $sp: 已清洗 ($clean_lines 行)" | tee -a $LOG
            skip=$((skip + 1))
            continue
        fi
    fi
    
    awk '/^>/{print $1; next} {print}' $raw > $clean.tmp
    if [ $? -ne 0 ]; then
        echo "  [$i/$n_species] ❌ $sp: 清洗失败" | tee -a $LOG
        \rm -f $clean.tmp
        fail=$((fail + 1))
        continue
    fi
    
    raw_lines=$(wc -l < $raw)
    tmp_lines=$(wc -l < $clean.tmp)
    if [ "$raw_lines" != "$tmp_lines" ]; then
        echo "  [$i/$n_species] ❌ $sp: 行数不符 raw=$raw_lines vs tmp=$tmp_lines" | tee -a $LOG
        \rm -f $clean.tmp
        fail=$((fail + 1))
        continue
    fi
    
    mv $clean.tmp $clean
    
    first_header=$(head -1 $clean)
    if echo "$first_header" | grep -q ' '; then
        echo "  [$i/$n_species] ❌ $sp: header 仍含空格: $first_header" | tee -a $LOG
        \rm -f $clean
        fail=$((fail + 1))
        continue
    fi
    
    sz=$(du -h $clean | awk '{print $1}')
    echo "  [$i/$n_species] ✅ $sp: $sz | $first_header" | tee -a $LOG
    ok=$((ok + 1))
done

echo "" | tee -a $LOG
echo "===== 清洗完成 =====" | tee -a $LOG
echo "  ✅ 成功清洗: $ok" | tee -a $LOG
echo "  ⏭  跳过（已存在）: $skip" | tee -a $LOG
echo "  ❌ 清洗失败: $fail" | tee -a $LOG
echo "  ❌ 原始缺失: $missing" | tee -a $LOG
echo "  合计: $((ok + skip + fail + missing)) / $n_species" | tee -a $LOG
echo "  结束时间: $(ts)" | tee -a $LOG
echo "" | tee -a $LOG

final_count=$(ls $GENOME_DIR/*_genomic_softmasked.clean.fa 2>/dev/null | wc -l)
echo "===== 最终 clean.fa 数量: $final_count =====" | tee -a $LOG
if [ "$final_count" -eq 125 ]; then
    echo "  ✅ 全部 125 株清洗完成，可以进入 Step 02 预检" | tee -a $LOG
    exit 0
else
    echo "  ⚠️  数量不符（期望 125），检查日志: $LOG" | tee -a $LOG
    exit 1
fi
