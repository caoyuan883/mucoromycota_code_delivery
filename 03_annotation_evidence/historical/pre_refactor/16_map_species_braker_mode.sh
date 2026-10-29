#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# 16_map_species_braker_mode.sh
# Builds a species-level mapping table that determines, for each softmasked
# genome, whether a merged RNA-seq BAM file is available (output of script
# 15), and assigns the corresponding BRAKER run mode accordingly
# (BRAKER3 if RNA-seq is available, BRAKER2 otherwise).
# ---------------------------------------------------------------------------

GENOME_DIR="${1:-/media/desk16/tly8529/mucoromycota_annotation_braker/03_softmasked_genomes}"
MERGED_DIR="${2:-/media/desk16/tly8529/mucoromycota_annotation_braker/RNA_SEQ/05_merged_bam}"
MAP_OUT="${3:-/media/desk16/tly8529/mucoromycota_annotation_braker/species_mode_mapping.csv}"

echo "species,rna_seq_status,merged_bam_file,braker_mode" > "$MAP_OUT"

for g in "$GENOME_DIR"/*_genomic_softmasked.fa; do
  base=$(basename "$g" _genomic_softmasked.fa)
  merged="$MERGED_DIR/${base}.merged.bam"
  if [ -f "$merged" ]; then
    echo "${base},Yes,${base}.merged.bam,BRAKER3" >> "$MAP_OUT"
  else
    echo "${base},No,NA,BRAKER2" >> "$MAP_OUT"
  fi
done

(head -n 1 "$MAP_OUT" && tail -n +2 "$MAP_OUT" | sort -t',' -k1,1) > "${MAP_OUT}.sorted"
mv "${MAP_OUT}.sorted" "$MAP_OUT"

echo "Mapping table written to: $MAP_OUT"
echo "----------------------------------------"
echo "Summary:"
echo "Species with RNA-seq (BRAKER3): $(grep -c 'BRAKER3' "$MAP_OUT")"
echo "Species without RNA-seq (BRAKER2): $(grep -c 'BRAKER2' "$MAP_OUT")"
echo "Total genomes: $(($(wc -l < "$MAP_OUT") - 1))"
