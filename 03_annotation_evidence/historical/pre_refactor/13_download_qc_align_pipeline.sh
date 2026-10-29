#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# 13_download_qc_align_pipeline.sh
# Downloads SRA runs, performs quality trimming (fastp), builds HISAT2
# genome indices, and aligns each run to its matched genome, producing
# sorted/indexed BAM files, alignment-rate summaries, and flagstat reports.
# ---------------------------------------------------------------------------

XLSX="${1:-18_tables_sra_run_download_list_with_genome_path.xlsx}"
BASE="${2:-rnaseq_braker_work}"

TABLES="$BASE/00_tables"
RAW="$BASE/01_fastq"
CLEAN="$BASE/02_clean"
IDXROOT="$BASE/03_index"
RUNBAM="$BASE/04_run_bam"
QC="$BASE/06_qc"
LOG="$BASE/07_logs"

mkdir -p "$TABLES" "$RAW" "$CLEAN" "$IDXROOT" "$RUNBAM" "$QC" "$LOG"

python3 - "$XLSX" > "$TABLES/run_rows.tsv" << 'PY'
import sys
import pandas as pd

xlsx_file = sys.argv[1]
df = pd.read_excel(xlsx_file)

for _, row in df.iterrows():
    species = str(row["Species"]).strip()
    run = str(row["Run_Accession"]).strip()
    layout = str(row["LibraryLayout"]).strip()
    genome = str(row["GenomePath"]).strip()
    species_tag = species.replace(" ", "_").replace("/", "_")
    print("\t".join([species, species_tag, run, layout, genome]))
PY

while IFS=$'\t' read -r SPECIES SPECIES_TAG RUN LAYOUT GENOME; do
    [ -z "$RUN" ] && continue

    RUN_BAM="$RUNBAM/${SPECIES_TAG}__${RUN}.sorted.bam"
    RUN_BAI="${RUN_BAM}.bai"
    RATE_FILE="$QC/${SPECIES_TAG}__${RUN}.alignment_rate.txt"
    STDERR_FILE="$QC/${SPECIES_TAG}__${RUN}.hisat2.stderr"
    FLAGSTAT_FILE="$QC/${SPECIES_TAG}__${RUN}.flagstat.txt"
    FASTP_JSON="$QC/${SPECIES_TAG}__${RUN}.fastp.json"
    FASTP_HTML="$QC/${SPECIES_TAG}__${RUN}.fastp.html"
    RUN_LOG="$LOG/${SPECIES_TAG}__${RUN}.log"

    IDXDIR="$IDXROOT/$SPECIES_TAG"
    IDXBASE="$IDXDIR/$SPECIES_TAG"

    mkdir -p "$IDXDIR"

    if [ -s "$RUN_BAM" ] && [ -s "$RUN_BAI" ] && [ -s "$RATE_FILE" ]; then
        echo "[SKIP] $SPECIES / $RUN already done"
        continue
    fi

    exec > >(tee -a "$RUN_LOG") 2>&1
    echo "[$(date)] Start $SPECIES / $RUN / $LAYOUT"

    if [ ! -s "${IDXBASE}.1.ht2" ]; then
        echo "[$(date)] Building HISAT2 index for $SPECIES"
        hisat2-build -p 4 "$GENOME" "$IDXBASE"
    fi

    cd "$RAW"
    echo "[$(date)] Downloading $RUN"
    fasterq-dump "$RUN" -O "$RAW" -e 4 -p

    if [ "$LAYOUT" = "PAIRED" ]; then
        echo "[$(date)] fastp $RUN (paired)"
        fastp \
          -i "$RAW/${RUN}_1.fastq" \
          -I "$RAW/${RUN}_2.fastq" \
          -o "$CLEAN/${SPECIES_TAG}__${RUN}_1.clean.fastq.gz" \
          -O "$CLEAN/${SPECIES_TAG}__${RUN}_2.clean.fastq.gz" \
          -w 4 \
          -h "$FASTP_HTML" \
          -j "$FASTP_JSON"

        echo "[$(date)] HISAT2 $RUN (paired)"
        hisat2 -p 4 -x "$IDXBASE" \
          -1 "$CLEAN/${SPECIES_TAG}__${RUN}_1.clean.fastq.gz" \
          -2 "$CLEAN/${SPECIES_TAG}__${RUN}_2.clean.fastq.gz" \
          2> "$STDERR_FILE" \
        | samtools sort -@ 2 -o "$RUN_BAM"

        rm -f "$RAW/${RUN}_1.fastq" "$RAW/${RUN}_2.fastq"
        rm -f "$CLEAN/${SPECIES_TAG}__${RUN}_1.clean.fastq.gz" "$CLEAN/${SPECIES_TAG}__${RUN}_2.clean.fastq.gz"
    else
        echo "[$(date)] fastp $RUN (single)"
        fastp \
          -i "$RAW/${RUN}.fastq" \
          -o "$CLEAN/${SPECIES_TAG}__${RUN}.clean.fastq.gz" \
          -w 4 \
          -h "$FASTP_HTML" \
          -j "$FASTP_JSON"

        echo "[$(date)] HISAT2 $RUN (single)"
        hisat2 -p 4 -x "$IDXBASE" \
          -U "$CLEAN/${SPECIES_TAG}__${RUN}.clean.fastq.gz" \
          2> "$STDERR_FILE" \
        | samtools sort -@ 2 -o "$RUN_BAM"

        rm -f "$RAW/${RUN}.fastq"
        rm -f "$CLEAN/${SPECIES_TAG}__${RUN}.clean.fastq.gz"
    fi

    samtools index "$RUN_BAM"
    samtools flagstat "$RUN_BAM" > "$FLAGSTAT_FILE"

    python3 - "$STDERR_FILE" "$RATE_FILE" << 'PY'
import re, sys
stderr, out = sys.argv[1], sys.argv[2]
rate = "NA"
with open(stderr, encoding="utf-8", errors="ignore") as f:
    text = f.read()
m = re.search(r'([0-9.]+)% overall alignment rate', text)
if m:
    rate = m.group(1)
with open(out, "w") as g:
    g.write(rate + "\n")
PY

    echo "[$(date)] Done $SPECIES / $RUN"
done < "$TABLES/run_rows.tsv"
