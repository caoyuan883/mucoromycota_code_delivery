#!/usr/bin/env bash
set -euo pipefail

B2="/home/py/Desktop/braker2_127_proteins"
B3="/home/py/Desktop/braker3_83_proteins"
OUT="/home/py/Desktop/Mucoromycota_BRAKER_core_proteomes"

[ -d "$B2" ] || { echo "找不到：$B2"; exit 1; }
[ -d "$B3" ] || { echo "找不到：$B3"; exit 1; }

N2=$(find "$B2" -maxdepth 1 -type f -name '*_protein.faa' | wc -l)
N3=$(find "$B3" -maxdepth 1 -type f -name '*_protein.faa' | wc -l)

echo "BRAKER2 文件数：$N2"
echo "BRAKER3 文件数：$N3"

DUP=$(comm -12 \
  <(find "$B2" -maxdepth 1 -type f -name '*_protein.faa' -printf '%f\n' | sort) \
  <(find "$B3" -maxdepth 1 -type f -name '*_protein.faa' -printf '%f\n' | sort))

if [ -n "$DUP" ]; then
    echo "发现同名文件，为防止覆盖，停止复制："
    echo "$DUP"
    exit 1
fi

mkdir -p "$OUT"

EXISTING=$(find "$OUT" -maxdepth 1 -type f -name '*_protein.faa' | wc -l)

if [ "$EXISTING" -gt 0 ]; then
    echo "输出目录中已经有 $EXISTING 个蛋白文件。"
    echo "为防止覆盖，本次停止。"
    exit 1
fi

find "$B2" -maxdepth 1 -type f -name '*_protein.faa' \
  -exec cp -p -t "$OUT" {} +

find "$B3" -maxdepth 1 -type f -name '*_protein.faa' \
  -exec cp -p -t "$OUT" {} +

NOUT=$(find "$OUT" -maxdepth 1 -type f -name '*_protein.faa' | wc -l)
EXPECTED=$((N2 + N3))

{
    printf "protein_file\tsource\n"

    find "$B2" -maxdepth 1 -type f -name '*_protein.faa' \
      -printf '%f\tBRAKER2\n'

    find "$B3" -maxdepth 1 -type f -name '*_protein.faa' \
      -printf '%f\tBRAKER3\n'
} | sort > "$OUT/protein_file_manifest.tsv"

echo "========================================"
echo "预期文件数：$EXPECTED"
echo "实际文件数：$NOUT"
echo "输出目录：$OUT"

if [ "$NOUT" -eq "$EXPECTED" ]; then
    echo "合并成功。"
else
    echo "文件数量不一致，请检查。"
    exit 1
fi
