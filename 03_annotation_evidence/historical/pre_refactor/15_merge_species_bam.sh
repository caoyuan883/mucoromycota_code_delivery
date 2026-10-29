#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# 15_merge_species_bam.sh
# Merges per-run sorted BAM files (output of script 13) into a single
# species-level BAM file, grouping runs by the species tag prefix in the
# filename (Species_tag__Run_Accession.sorted.bam). Species with a single
# run are simply copied through; species with multiple runs are merged,
# re-sorted, and indexed.
# ---------------------------------------------------------------------------

BAM_DIR="${1:-/media/desk16/tly8529/mucoromycota_annotation_braker/RNA_SEQ/04_run_bam}"
MERGED_DIR="${2:-/media/desk16/tly8529/mucoromycota_annotation_braker/RNA_SEQ/05_merged_bam}"

mkdir -p "$MERGED_DIR"

unset species_bams
declare -A species_bams

cd "$BAM_DIR"
for f in *.sorted.bam; do
  sp="${f%%__*}"
  species_bams["$sp"]+="$f "
done

for sp in $(printf '%s\n' "${!species_bams[@]}" | sort); do
  bams=($(printf '%s\n' ${species_bams[$sp]} | sort -u))
  out="$MERGED_DIR/${sp}.merged.bam"

  if [ "${#bams[@]}" -eq 1 ]; then
    cp "${bams[0]}" "$out"
    echo "[SINGLE] $sp -> ${bams[0]}"
  else
    samtools merge -f -@ 8 -o "${out}.tmp.bam" "${bams[@]}"
    samtools sort -@ 8 -o "$out" "${out}.tmp.bam"
    rm "${out}.tmp.bam"
    echo "[MERGED] $sp -> ${bams[*]}"
  fi
  samtools index -@ 4 "$out"
done

echo "----------------------------------------"
echo "共处理物种数: ${#species_bams[@]}"
