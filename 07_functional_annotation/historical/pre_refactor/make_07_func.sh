#!/bin/bash
set -uo pipefail
ROOT=/Users/caoyuan/Desktop/interproscan_eggnog
EMIN=$ROOT/func_freeze_2/emapper_v1_20261011
IPSIN=$ROOT/func_freeze_2/interproscan_v1_20261011
EM10=$ROOT/func_filtered10_release/emapper_v1_20261011_filtered10
IPS10=$ROOT/func_filtered10_release/interproscan_v1_20261011_filtered10
FAADIR=/Users/caoyuan/Desktop/Mucoromycota_BRAKER_longest_isoform/02_202_longest_isoform_proteins
OUT=$ROOT/07_202_functional_annotation
[ -e "$OUT" ] && { echo "已存在 $OUT，请先改名或删除"; exit 1; }
mkdir -p "$OUT"
LOG=$OUT/00_build_log.tsv
TMP=$(mktemp -d)
printf "strain\tsource\tn_faa\tegg_proteins\tegg_not_in_faa\tips_rows\tips_proteins\tips_not_in_faa\n" > "$LOG"

for faa in "$FAADIR"/*.faa; do
  sp=$(basename "$faa"); sp=${sp%%_longest_isoform*}
  if [ -f "$EM10/${sp}.emapper.annotations" ]; then
    src=filtered10; em=$EM10/${sp}.emapper.annotations; ips=$IPS10/${sp}.interproscan.tsv.gz
  else
    src=freeze2;    em=$EMIN/${sp}.emapper.annotations;  ips=$IPSIN/${sp}.interproscan.tsv.gz
  fi
  if [ ! -f "$em" ] || [ ! -f "$ips" ]; then echo "缺文件: $sp ($src)"; continue; fi
  d=$OUT/$sp; mkdir -p "$d"
  grep -v '^##' "$em" > "$d/${sp}_emapper.annotations.tsv"
  gzcat "$ips" > "$d/${sp}_interproscan.tsv"
  grep '>' "$faa" | sed 's/^>//;s/[[:space:]].*//' | sort -u > "$TMP/ids"
  awk -F'\t' '!/^#/{print $1}' "$d/${sp}_emapper.annotations.tsv" | sort -u > "$TMP/em"
  cut -f1 "$d/${sp}_interproscan.tsv" | sort -u > "$TMP/ips"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sp" "$src" \
    "$(wc -l < "$TMP/ids" | tr -d ' ')" \
    "$(wc -l < "$TMP/em" | tr -d ' ')"  "$(comm -23 "$TMP/em"  "$TMP/ids" | wc -l | tr -d ' ')" \
    "$(wc -l < "$d/${sp}_interproscan.tsv" | tr -d ' ')" \
    "$(wc -l < "$TMP/ips" | tr -d ' ')" "$(comm -23 "$TMP/ips" "$TMP/ids" | wc -l | tr -d ' ')" >> "$LOG"
done
rm -rf "$TMP"

echo "== 汇总 =="
awk -F'\t' 'NR>1{n++; if($2=="filtered10")f++; a+=$3; e+=$4; x+=$5+$8; i+=$7}
  END{printf "株数=%d (期望202)  过滤株=%d (期望10)\n蛋白=%d (期望2618164)\negg有注释=%d (期望2143996)\nips有命中=%d (期望2000186)\n注释ID不在faa中=%d (期望0)\n",n,f,a,e,i,x}' "$LOG"

cd "$OUT" && find . -mindepth 2 -type f -name "*.tsv" | sort | xargs md5 -r > MD5SUMS.txt
echo "完成：$OUT"
