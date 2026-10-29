#!/bin/bash
set -uo pipefail
B=/Users/caoyuan/Desktop/Mucoromycota_BRAKER_longest_isoform
GTF=$B/01_202_longest_isoform_gtf
FAA=$B/02_202_longest_isoform_proteins
GFF=$B/06_202_longest_isoform_gff3
FUN=/Users/caoyuan/Desktop/interproscan_eggnog/07_202_functional_annotation
OUT=/Users/caoyuan/Desktop/Mucoromycota_202_release
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
awk -F'\t' 'NR>1{n++; p+=$2; d+=$5+$6+$7; if($8!="none")m++}
  END{printf "株数=%d (期望202)\n蛋白总数=%d (期望2618164)\nID不一致合计=%d (期望0)\n缺文件的株=%d (期望0)\n",n,p,d,m}' "$LOG"
echo "每株文件数不是5的："
for d in "$OUT"/*/; do c=$(ls "$d" | wc -l | tr -d ' '); [ "$c" -ne 5 ] && echo "  $(basename "$d") $c"; done
cd "$OUT" && find . -mindepth 2 -type f | sort | xargs md5 -r > MD5SUMS.txt
echo "完成：$OUT"
