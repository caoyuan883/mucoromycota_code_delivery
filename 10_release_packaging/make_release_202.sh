#!/bin/bash
# Assemble the per-strain release folders (FAA, GTF, GFF3, InterProScan TSV, eggNOG TSV).
# [release refactor] Paths come from the environment; copy and check logic unchanged.
#   REP_ROOT=/path/to/Mucoromycota_BRAKER_longest_isoform \
#   FUNC_RELEASE=/path/to/07_202_functional_annotation \
#   OUT=/path/to/Mucoromycota_202_release  bash make_release_202.sh
#   STRICT=1 ... -> exit 3 if the printed expectations are not met
# Per-file checksums are written to MD5SUMS_files.txt (historically MD5SUMS.txt, which
# pack_release_202.sh later overwrote with the archive checksums in the same folder).
set -uo pipefail
# historical: B=/Users/caoyuan/Desktop/Mucoromycota_BRAKER_longest_isoform
B=${REP_ROOT:?set REP_ROOT (folder with 01_202_longest_isoform_gtf, 02_202_longest_isoform_proteins, 06_202_longest_isoform_gff3)}
GTF=${GTF:-$B/01_202_longest_isoform_gtf}
FAA=${FAA:-$B/02_202_longest_isoform_proteins}
GFF=${GFF:-$B/06_202_longest_isoform_gff3}
# historical: FUN=/Users/caoyuan/Desktop/interproscan_eggnog/07_202_functional_annotation
FUN=${FUNC_RELEASE:?set FUNC_RELEASE (output folder of assemble_functional_tables.sh)}
# historical: OUT=/Users/caoyuan/Desktop/Mucoromycota_202_release
OUT=${OUT:?set OUT (new release folder; must not exist)}
STRICT=${STRICT:-0}
for d in "$GTF" "$FAA" "$GFF" "$FUN"; do
  [ -d "$d" ] || { echo "ERROR: directory not found: $d" >&2; exit 2; }
done
md5list(){ if command -v md5sum >/dev/null 2>&1; then md5sum "$@"; else md5 -r "$@"; fi; }
[ -e "$OUT" ] && { echo "已存在 $OUT，请先改名或删除"; exit 1; }
mkdir -p "$OUT"
LOG=$OUT/00_release_check.tsv
TMP=$(mktemp -d)
printf "strain\tn_faa\tn_gtf_tx\tn_gff3_mrna\tfaa_not_in_gtf\tgtf_not_in_faa\tfaa_not_in_gff3\tmissing\n" > "$LOG"

for faa in "$FAA"/*.faa; do
  sp=$(basename "$faa"); sp=${sp%%_longest_isoform*}
  gtf=$GTF/${sp}_longest_isoform.gtf
  gff=$(find "$GFF" -type f -name "${sp}_longest_isoform.gff3" | head -1)
  ips=$FUN/$sp/${sp}_interproscan.tsv
  em=$FUN/$sp/${sp}_emapper.annotations.tsv
  miss=""
  for f in "$gtf" "$gff" "$ips" "$em"; do [ -s "$f" ] || miss="$miss $(basename "${f:-gff3}")"; done
  d=$OUT/$sp; mkdir -p "$d"
  cp -p "$faa" "$d/${sp}_longest_isoform.faa"
  [ -s "$gtf" ] && cp -p "$gtf" "$d/${sp}_longest_isoform.gtf"
  [ -s "$gff" ] && cp -p "$gff" "$d/${sp}_longest_isoform.gff3"
  [ -s "$ips" ] && cp -p "$ips" "$d/${sp}_interproscan.tsv"
  [ -s "$em" ]  && cp -p "$em"  "$d/${sp}_emapper.annotations.tsv"

  grep '>' "$faa" | sed 's/^>//;s/[[:space:]].*//' | sort -u > "$TMP/faa"
  awk -F'\t' 'match($9,/transcript_id "[^"]+"/){print substr($9,RSTART+15,RLENGTH-16)}' "$gtf" 2>/dev/null | sort -u > "$TMP/gtf"
  awk -F'\t' '($3=="mRNA"||$3=="transcript") && match($9,/ID=[^;]+/){print substr($9,RSTART+3,RLENGTH-3)}' "$gff" 2>/dev/null | sort -u > "$TMP/gff"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sp" \
    "$(wc -l < "$TMP/faa" | tr -d ' ')" "$(wc -l < "$TMP/gtf" | tr -d ' ')" "$(wc -l < "$TMP/gff" | tr -d ' ')" \
    "$(comm -23 "$TMP/faa" "$TMP/gtf" | wc -l | tr -d ' ')" "$(comm -13 "$TMP/faa" "$TMP/gtf" | wc -l | tr -d ' ')" \
    "$(comm -23 "$TMP/faa" "$TMP/gff" | wc -l | tr -d ' ')" "${miss:-none}" >> "$LOG"
done
rm -rf "$TMP"

echo "== 汇总 =="
awk -F'\t' -v strict="$STRICT" 'NR>1{n++; p+=$2; d+=$5+$6+$7; if($8!="none")m++}
  END{printf "株数=%d (期望202)\n蛋白总数=%d (期望2618164)\nID不一致合计=%d (期望0)\n缺文件的株=%d (期望0)\n",n,p,d,m
      if (strict==1 && (n!=202 || p!=2618164 || d!=0 || m!=0)) { print "STRICT: expectation mismatch" > "/dev/stderr"; exit 3 } }' "$LOG"
rc=$?
echo "每株文件数不是5的："
bad5=0
for d in "$OUT"/*/; do c=$(ls "$d" | wc -l | tr -d ' '); [ "$c" -ne 5 ] && { echo "  $(basename "$d") $c"; bad5=$((bad5+1)); }; done
cd "$OUT" && find . -mindepth 2 -type f | sort | while IFS= read -r f; do md5list "$f"; done > MD5SUMS_files.txt
echo "完成：$OUT"
if [ "$STRICT" = 1 ] && { [ "$rc" -ne 0 ] || [ "$bad5" -ne 0 ]; }; then exit 3; fi
