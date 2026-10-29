#!/bin/bash
# ============================================================================
# Step 0-a: 从 species_mode_mapping.csv 提取 86 株 BRAKER3 物种清单
# 输出: species_list_braker3.txt (每行一个物种名)
# ============================================================================
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)

CSV=/media/desk16/tly8529/mucoromycota_annotation_braker/species_mode_mapping.csv
OUT=$here/species_list_braker3.txt

if [[ ! -s "$CSV" ]]; then
  echo "❌ CSV 不存在: $CSV" >&2
  exit 1
fi

awk -F',' '$4=="BRAKER3" {print $1}' "$CSV" > "$OUT"
N=$(wc -l < "$OUT")

echo "✅ 已生成 $OUT"
echo "   总数: $N 株"
echo "   前 5 株: $(head -5 "$OUT" | tr '\n' ' ')"
echo "   后 5 株: $(tail -5 "$OUT" | tr '\n' ' ')"

if [[ "$N" -ne 86 ]]; then
  echo "⚠️  数量不是 86，请确认 CSV 内容是否正确" >&2
fi
