#!/bin/bash
# ============================================================================
# Step 0-b: 批量 header 清洗
# 规则: 保留 ">" 后第一个空格前的字符串（GenBank accession）
# 例子:
#   >MCGE01000001.1 Absidia repens strain NRRL 1336 ...  →  >MCGE01000001.1
#
# 输入: $GENOME_DIR/<sp>_genomic_softmasked.fa
# 输出: $GENOME_DIR/<sp>_genomic_softmasked.clean.fa
#
# 安全策略：
#   - 已存在 .clean.fa 且非空 → 跳过（不覆盖）
#   - 用 awk 一次性生成，绝不 in-place 修改原文件
#   - 全部完成后打印统计
# ============================================================================
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

LIST=$here/species_list_braker3.txt
if [[ ! -s "$LIST" ]]; then
  echo "❌ $LIST 不存在，请先跑 prepare_species_list.sh" >&2
  exit 1
fi

CLEAN_LOG=$here/00b_clean_headers.log
: > "$CLEAN_LOG"

n_total=0; n_skip=0; n_clean=0; n_fail=0
while IFS= read -r SP; do
  [[ -z "$SP" ]] && continue
  n_total=$((n_total+1))
  RAW=$GENOME_DIR/${SP}_genomic_softmasked.fa
  CLEAN=$GENOME_DIR/${SP}_genomic_softmasked.clean.fa

  if [[ ! -s "$RAW" ]]; then
    echo "[$SP] ❌ 原始 fa 不存在: $RAW" | tee -a "$CLEAN_LOG"
    n_fail=$((n_fail+1))
    continue
  fi

  if [[ -s "$CLEAN" ]]; then
    # 已存在，且大小合理
    n_skip=$((n_skip+1))
    echo "[$SP] ⏭  已存在 clean.fa，跳过" | tee -a "$CLEAN_LOG"
    continue
  fi

  # awk 生成：header 只保留第一个字段（空格分割）
  # $1 已经包含 ">"，直接输出
  if awk '/^>/{print $1; next} {print}' "$RAW" > "$CLEAN.tmp"; then
    mv "$CLEAN.tmp" "$CLEAN"
    n_headers_before=$(grep -c '^>' "$RAW")
    n_headers_after=$(grep -c '^>' "$CLEAN")
    # 校验：header 数应一致
    if [[ "$n_headers_before" -eq "$n_headers_after" ]]; then
      # 再校验：clean.fa 里不应还有空格 header
      n_bad=$(grep '^>' "$CLEAN" | grep -c ' ' || true)
      if [[ "$n_bad" -eq 0 ]]; then
        n_clean=$((n_clean+1))
        echo "[$SP] ✅ 已清洗 (n_seq=$n_headers_after)" | tee -a "$CLEAN_LOG"
      else
        echo "[$SP] ❌ 清洗后仍有 $n_bad 条含空格 header" | tee -a "$CLEAN_LOG"
        n_fail=$((n_fail+1))
      fi
    else
      echo "[$SP] ❌ header 数不一致 before=$n_headers_before after=$n_headers_after" | tee -a "$CLEAN_LOG"
      n_fail=$((n_fail+1))
    fi
  else
    echo "[$SP] ❌ awk 处理失败" | tee -a "$CLEAN_LOG"
    n_fail=$((n_fail+1))
    /bin/rm -f "$CLEAN.tmp"
  fi
done < "$LIST"

echo ""
echo "========================================================"
echo "Header 清洗汇总"
echo "========================================================"
echo "总数:     $n_total"
echo "新清洗:   $n_clean"
echo "已存在:   $n_skip"
echo "失败:     $n_fail"
echo "日志:     $CLEAN_LOG"
echo "========================================================"

if [[ "$n_fail" -gt 0 ]]; then
  exit 2
fi
