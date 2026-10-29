#!/usr/bin/env bash
set -uo pipefail

# [release refactor] historical: BASE=/home/py/Desktop/12_faa_12_gff_technical_validation
BASE="${REF12_BASE:?set REF12_BASE to the 12-genome technical-validation folder}"
REF="$BASE/12_faa_longest_isoform"
OWN="$BASE/12_faa_own_longest_isoform"
WORK="$BASE/protein_comparison"
mkdir -p "$WORK/db" "$WORK/aln" "$WORK/log" "$WORK/tmp"

command -v diamond >/dev/null || { echo "错误：找不到 diamond"; exit 1; }
echo "DIAMOND: $(diamond --version 2>&1 | head -1)"

# 线程数自动取，最多 8，避免小机器超订
NP=$(nproc); T=$(( NP > 8 ? 8 : NP )); [ "$T" -lt 1 ] && T=1
# 低内存模式：约 2.5 GB 峰值
BLK="0.4"; CHK="1"
FMT="qseqid sseqid pident length qlen slen qstart qend sstart send bitscore evalue"
echo "线程 $T   block $BLK   chunks $CHK"
echo

# 只跑某一株：  bash run_protein_comparison.sh Lentamyces_parricidus
ONLY="${1:-}"

run_one () {
  local q="$1" db="$2" out="$3" lab="$4"
  if [ -s "$out" ]; then
    echo "    SKIP $lab (已有 $(wc -l < "$out") 行)"; return 0
  fi
  echo "    [$(date +%H:%M:%S)] $lab  开始"
  diamond blastp \
    -q "$q" -d "$db" -o "$out.partial" \
    --outfmt 6 $FMT \
    --max-target-seqs 10 --evalue 1e-5 \
    -p "$T" -b "$BLK" -c "$CHK" \
    --tmpdir "$WORK/tmp" \
    --log 2>&1 | tee -a "$WORK/log/${lab}.log" | grep -E 'Processing|Total time|Reported|Error' 
  local rc=${PIPESTATUS[0]}
  if [ "$rc" -eq 0 ]; then
    mv "$out.partial" "$out"
    echo "    [$(date +%H:%M:%S)] $lab  完成  $(wc -l < "$out") 行"
  else
    rm -f "$out.partial"
    echo "    [$(date +%H:%M:%S)] $lab  失败 rc=$rc  见 $WORK/log/${lab}.log"
    return 1
  fi
}

# 按文件大小从小到大排，先跑小的
for f in $(ls -S -r "$REF"/*_longest_isoform.faa); do
  s=$(basename "$f" _longest_isoform.faa)
  [ -n "$ONLY" ] && [ "$s" != "$ONLY" ] && continue
  o="$OWN/${s}_longest_isoform.faa"
  [ -f "$o" ] || { echo "[缺本方] $s"; continue; }

  echo "=== $s   ref $(grep -c '^>' "$f") 条 / own $(grep -c '^>' "$o") 条 ==="

  for side in ref own; do
    src=$([ "$side" = ref ] && echo "$f" || echo "$o")
    if [ ! -f "$WORK/db/${s}_${side}.dmnd" ]; then
      echo "    建库 $side"
      diamond makedb --in "$src" -d "$WORK/db/${s}_${side}" -p "$T" \
        >> "$WORK/log/${s}.makedb.log" 2>&1 \
        || { echo "    建库失败，见 $WORK/log/${s}.makedb.log"; continue 2; }
    fi
  done

  run_one "$o" "$WORK/db/${s}_ref" "$WORK/aln/${s}.own_vs_ref.tsv" "${s}.own_vs_ref"
  run_one "$f" "$WORK/db/${s}_own" "$WORK/aln/${s}.ref_vs_own.tsv" "${s}.ref_vs_own"
  echo
done
echo "[$(date +%H:%M:%S)] 全部结束"
