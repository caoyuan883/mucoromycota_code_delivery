#!/bin/bash
# Convert representative-isoform GTF files to structural GFF3 (06_202_longest_isoform_gff3).
#
# Status: reconstructed. The original batch wrapper was not preserved; the per-file
# command used for the release is:
#     perl "$G2G" --gff3 --out="$output_gff3" < "$input_gtf"
# This wrapper runs exactly that command, without --printExon/--printIntron
# (unlike the BRAKER2 raw-output helper 04_structural_annotation/braker2/make_raw_gff3.sh).
# The historical converter revision was not recorded; record the one you use.
#
# Usage:
#   G2G=/path/to/Augustus/scripts/gtf2gff.pl \
#   GTF_DIR=/path/to/01_202_longest_isoform_gtf \
#   GFF3_DIR=/path/to/new_gff3_folder \
#   bash gtf_to_gff3.sh
# Optional: GTF_SUFFIX (default _longest_isoform.gtf), GFF3_SUFFIX (default _longest_isoform.gff3)
# Exit codes: 0 all converted; 1 output folder not empty; 2 missing input; 3 conversion failure(s)
set -u
G2G=${G2G:?set G2G (path to gtf2gff.pl)}
GTF_DIR=${GTF_DIR:?set GTF_DIR (folder with representative GTF files)}
GFF3_DIR=${GFF3_DIR:?set GFF3_DIR (new output folder)}
GTF_SUFFIX=${GTF_SUFFIX:-_longest_isoform.gtf}
GFF3_SUFFIX=${GFF3_SUFFIX:-_longest_isoform.gff3}

[ -f "$G2G" ] || { echo "ERROR: gtf2gff.pl not found: $G2G" >&2; exit 2; }
[ -d "$GTF_DIR" ] || { echo "ERROR: GTF folder not found: $GTF_DIR" >&2; exit 2; }
if [ -d "$GFF3_DIR" ] && [ -n "$(ls -A "$GFF3_DIR" 2>/dev/null)" ]; then
  echo "ERROR: output folder not empty: $GFF3_DIR" >&2; exit 1
fi
mkdir -p "$GFF3_DIR"
LOG=$GFF3_DIR/00_gtf_to_gff3_log.tsv
printf "label\tinput_gtf\toutput_gff3\texit_code\tstatus\n" > "$LOG"
{ echo "gtf2gff.pl: $G2G"; sha256sum "$G2G" 2>/dev/null || shasum -a 256 "$G2G"; perl -v | sed -n 2p; } > "$GFF3_DIR/00_tool_version.txt"

n=0; fail=0
shopt -s nullglob
for input_gtf in "$GTF_DIR"/*"$GTF_SUFFIX"; do
  label=$(basename "$input_gtf" "$GTF_SUFFIX")
  output_gff3=$GFF3_DIR/${label}${GFF3_SUFFIX}
  perl "$G2G" --gff3 --out="$output_gff3" < "$input_gtf"
  rc=$?
  if [ $rc -eq 0 ] && [ -s "$output_gff3" ]; then st=OK; else st=FAIL; fail=$((fail+1)); fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$label" "$input_gtf" "$output_gff3" "$rc" "$st" >> "$LOG"
  n=$((n+1))
done
[ $n -gt 0 ] || { echo "ERROR: no *$GTF_SUFFIX files in $GTF_DIR" >&2; exit 2; }
echo "converted=$((n-fail)) failed=$fail  log: $LOG"
[ $fail -eq 0 ] || exit 3
