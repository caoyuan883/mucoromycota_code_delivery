#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# summarize_alignment_qc.sh
# Aggregates per-run alignment rate, flagstat, and fastp QC metrics
# (outputs of script 13) into a single summary table.
# ---------------------------------------------------------------------------

# [release refactor] historical default: /media/desk16/tly8529/mucoromycota_annotation_braker/RNA_SEQ
BASE="${1:?usage: $0 RNA_SEQ_DIR}"
BAM_DIR="$BASE/04_run_bam"
QC_DIR="$BASE/06_qc"
OUT="$BASE/rnaseq_qc_summary.tsv"

echo -e "Sample\tAlignmentRate(%)\tTotalReads\tMappedReads\tProperlyPaired\tDuplicationRate(%)\tReadsAfterFilter\tInsertSizePeak\tBAM_Size_MB" > "$OUT"

for RATE_FILE in "$QC_DIR"/*.alignment_rate.txt; do
    SAMPLE=$(basename "$RATE_FILE" .alignment_rate.txt)
    RATE=$(cat "$RATE_FILE" 2>/dev/null | tr -d '\n')

    FLAGSTAT="$QC_DIR/${SAMPLE}.flagstat.txt"
    TOTAL=$(grep -m1 "in total" "$FLAGSTAT" 2>/dev/null | awk '{print $1}')
    MAPPED=$(grep -m1 " mapped (" "$FLAGSTAT" 2>/dev/null | awk '{print $1}')
    PROPER=$(grep -m1 "properly paired" "$FLAGSTAT" 2>/dev/null | awk '{print $1}')

    JSON="$QC_DIR/${SAMPLE}.fastp.json"
    if [ -f "$JSON" ]; then
        DUP=$(grep -o '"duplication_rate": [0-9.]*' "$JSON" | head -1 | awk -F': ' '{print $2}')
        AFTER=$(grep -A2 '"after_filtering"' "$JSON" | grep '"total_reads"' | head -1 | grep -o '[0-9]*')
        PEAK=$(grep -o '"peak": [0-9]*' "$JSON" | head -1 | awk -F': ' '{print $2}')
    else
        DUP="NA"; AFTER="NA"; PEAK="NA"
    fi

    BAM="$BAM_DIR/${SAMPLE}.sorted.bam"
    if [ -f "$BAM" ]; then
        SIZE=$(du -m "$BAM" | cut -f1)
    else
        SIZE="NA"
    fi

    echo -e "${SAMPLE}\t${RATE:-NA}\t${TOTAL:-NA}\t${MAPPED:-NA}\t${PROPER:-NA}\t${DUP:-NA}\t${AFTER:-NA}\t${PEAK:-NA}\t${SIZE:-NA}" >> "$OUT"
done

echo "Summary table written to: $OUT"
