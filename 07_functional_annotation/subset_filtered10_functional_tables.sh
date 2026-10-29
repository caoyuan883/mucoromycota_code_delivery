#!/bin/bash
# Subset frozen eggNOG/InterProScan outputs to the protein IDs kept in the 10 post-filter FAA files.
# [release refactor] Paths come from the environment (historical values in comments); logic unchanged.
#   FUNC_ROOT=/path/to/interproscan_eggnog FILTDIR=/path/to/10_omark_filtered_faa bash subset_filtered10_functional_tables.sh
#   STRICT=1 ... -> exit 3 if any DIFF against expected_from_02.tsv
# Note: the eggNOG file search keeps the historical pattern "*${sp}.emapper.*" (leading wildcard);
# among the 202 release labels no label is a suffix of another (checked for this release).
set -uo pipefail
# historical: ROOT=/Users/caoyuan/Desktop/interproscan_eggnog
ROOT=${FUNC_ROOT:?set FUNC_ROOT (directory containing func_freeze_2/)}
EMIN=${EMIN:-$ROOT/func_freeze_2/emapper_v1_20261011}
IPSIN=${IPSIN:-$ROOT/func_freeze_2/interproscan_v1_20261011}
# historical: FILTDIR=/Users/caoyuan/Desktop/10_omark_filtered_faa
FILTDIR=${FILTDIR:?set FILTDIR (directory with <strain>_longest_isoform.filtered.faa)}
OUT=${OUT:-$ROOT/func_filtered10_release}
STRICT=${STRICT:-0}
for d in "$EMIN" "$IPSIN" "$FILTDIR"; do
  [ -d "$d" ] || { echo "ERROR: directory not found: $d" >&2; exit 2; }
done
md5list(){ if command -v md5sum >/dev/null 2>&1; then md5sum "$@"; else md5 -r "$@"; fi; }
EMOUT=$OUT/emapper_v1_20261011_filtered10
IPSOUT=$OUT/interproscan_v1_20261011_filtered10
[ -e "$OUT" ] && { echo "已存在 $OUT，请先改名或删除"; exit 1; }
mkdir -p "$EMOUT" "$IPSOUT" "$OUT/keep_ids"
LOG=$OUT/00_filter_log.tsv
printf "strain\tfile\tremoved_ids_in_input\tlines_in\tlines_out\tproteins_out\n" > "$LOG"

rd(){ if [[ $1 == *.gz ]]; then gzip -dc "$1"; else cat "$1"; fi; }

for FILT in "$FILTDIR"/*_longest_isoform.filtered.faa; do
  sp=$(basename "$FILT" _longest_isoform.filtered.faa)
  keep=$OUT/keep_ids/${sp}.keep.txt
  grep '>' "$FILT" | sed 's/^>//;s/[[:space:]].*//' | sort -u > "$keep"
  n=0
  while IFS= read -r f; do
    [ -z "$f" ] && continue
    n=$((n+1))
    case "$f" in "$EMIN"/*) od=$EMOUT;; *) od=$IPSOUT;; esac
    o=$od/$(basename "$f")
    rem=$(rd "$f" | awk -F'\t' '!/^#/{print $1}' | sort -u | comm -23 - "$keep" | wc -l | tr -d ' ')
    if [[ $f == *.gz ]]; then
      rd "$f" | awk -F'\t' 'NR==FNR{k[$1];next} /^#/||($1 in k)' "$keep" - | gzip -n > "$o"
    else
      rd "$f" | awk -F'\t' 'NR==FNR{k[$1];next} /^#/||($1 in k)' "$keep" - > "$o"
    fi
    li=$(rd "$f" | awk '!/^#/' | wc -l | tr -d ' ')
    lo=$(rd "$o" | awk '!/^#/' | wc -l | tr -d ' ')
    po=$(rd "$o" | awk -F'\t' '!/^#/{print $1}' | sort -u | wc -l | tr -d ' ')
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sp" "$(basename "$f")" "$rem" "$li" "$lo" "$po" >> "$LOG"
  done < <(find "$EMIN" -maxdepth 1 -type f -name "*${sp}.emapper.*"; \
           find "$IPSIN" -maxdepth 1 -type f -name "${sp}.interproscan.tsv.gz")
  echo "$sp: 处理 $n 个文件"
done

cat > "$OUT/expected_from_02.tsv" <<'T'
strain	egg_hit	ips_hit
Acaulospora_colombiana	8057	7703
Apophysomyces_sp._BC1015	10895	9672
Apophysomyces_sp._BC1034	10636	9426
Gigaspora_margarita	13347	12717
Gigaspora_rosea	14137	13207
Gryganskiella_cystojenkinii	10231	10423
Mortierella_sp._GBA39	12197	12075
Oehlia_diaphana	12446	11170
Podila_clonocystis	10598	10335
Rhizophagus_prolifer	19823	16514
T

echo; echo "== 与 02 最终表核对 =="
awk -F'\t' -v strict="$STRICT" 'NR==FNR{if(FNR>1){e[$1]=$2;i[$1]=$3};next}
  FNR>1{ if($2~/emapper\.annotations/){x=e[$1];t="egg"} else if($2~/interproscan/){x=i[$1];t="ips"} else next;
         printf "%-30s %s 得到=%s 02表=%s %s\n",$1,t,$6,x,($6==x?"PASS":"DIFF"); if($6!=x) bad++ }
  END{ if (strict==1 && bad>0) exit 3 }' \
  "$OUT/expected_from_02.tsv" "$LOG"
rc=$?

cd "$OUT" && find . -type f ! -name MD5SUMS.txt | sort | while IFS= read -r f; do md5list "$f"; done > MD5SUMS.txt
echo; echo "完成：$OUT"
if [ "$STRICT" = 1 ] && [ "$rc" -ne 0 ]; then echo "STRICT: DIFF against expected_from_02.tsv" >&2; exit 3; fi
