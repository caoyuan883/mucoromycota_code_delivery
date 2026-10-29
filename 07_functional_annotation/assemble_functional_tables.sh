#!/bin/bash
# Assemble the 202 per-strain release functional tables (eggNOG + InterProScan).
# [release refactor] Paths come from the environment (historical values in comments);
# selection/copy logic unchanged. Usage:
#   FUNC_ROOT=/path/to/interproscan_eggnog FAADIR=/path/to/release_faa bash assemble_functional_tables.sh
#   STRICT=1 ...   -> exit 3 when a strain lacks a result file or a printed expectation differs
set -uo pipefail
# historical: ROOT=/Users/caoyuan/Desktop/interproscan_eggnog
ROOT=${FUNC_ROOT:?set FUNC_ROOT (directory containing func_freeze_2/ and func_filtered10_release/)}
EMIN=${EMIN:-$ROOT/func_freeze_2/emapper_v1_20261011}
IPSIN=${IPSIN:-$ROOT/func_freeze_2/interproscan_v1_20261011}
EM10=${EM10:-$ROOT/func_filtered10_release/emapper_v1_20261011_filtered10}
IPS10=${IPS10:-$ROOT/func_filtered10_release/interproscan_v1_20261011_filtered10}
# historical: FAADIR=/Users/caoyuan/Desktop/Mucoromycota_BRAKER_longest_isoform/02_202_longest_isoform_proteins
# (must hold the FINAL release proteins: 192 unchanged + 10 post-filter)
FAADIR=${FAADIR:?set FAADIR (final release <strain>_longest_isoform*.faa directory)}
OUT=${OUT:-$ROOT/07_202_functional_annotation}
STRICT=${STRICT:-0}
for d in "$EMIN" "$IPSIN" "$EM10" "$IPS10" "$FAADIR"; do
  [ -d "$d" ] || { echo "ERROR: directory not found: $d" >&2; exit 2; }
done
# portable helpers (historical script used macOS gzcat and md5 -r)
zc(){ gzip -dc "$@"; }
md5list(){ if command -v md5sum >/dev/null 2>&1; then md5sum "$@"; else md5 -r "$@"; fi; }
nmiss=0
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
  if [ ! -f "$em" ] || [ ! -f "$ips" ]; then echo "缺文件: $sp ($src)"; nmiss=$((nmiss+1)); continue; fi
  d=$OUT/$sp; mkdir -p "$d"
  grep -v '^##' "$em" > "$d/${sp}_emapper.annotations.tsv"
  zc "$ips" > "$d/${sp}_interproscan.tsv"
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
awk -F'\t' -v strict="$STRICT" 'NR>1{n++; if($2=="filtered10")f++; a+=$3; e+=$4; x+=$5+$8; i+=$7}
  END{printf "株数=%d (期望202)  过滤株=%d (期望10)\n蛋白=%d (期望2618164)\negg有注释=%d (期望2143996)\nips有命中=%d (期望2000186)\n注释ID不在faa中=%d (期望0)\n",n,f,a,e,i,x
      if (strict==1 && (n!=202 || f!=10 || a!=2618164 || e!=2143996 || i!=2000186 || x!=0)) { print "STRICT: expectation mismatch" > "/dev/stderr"; exit 3 } }' "$LOG"
rc=$?

cd "$OUT" && find . -mindepth 2 -type f -name "*.tsv" | sort | while IFS= read -r f; do md5list "$f"; done > MD5SUMS.txt
echo "完成：$OUT"
if [ "$STRICT" = 1 ]; then
  [ "$nmiss" -eq 0 ] || { echo "STRICT: $nmiss strain(s) missing result files" >&2; exit 3; }
  [ "$rc" -eq 0 ] || exit 3
fi
