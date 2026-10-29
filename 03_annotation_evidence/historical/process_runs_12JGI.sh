#!/usr/bin/env bash
set -euo pipefail

CSV="$HOME/Desktop/rnaseq_braker_work_12_JGI/00_tables/12JGI_sci_grade_download_list_v2_with_genome.csv"
BASE="$HOME/Desktop/rnaseq_braker_work_12_JGI"

RAW="$BASE/01_fastq"
CLEAN="$BASE/02_clean"
IDXROOT="$BASE/03_index"
RUNBAM="$BASE/04_run_bam"
QC="$BASE/06_qc"
LOG="$BASE/07_logs"

mkdir -p "$RAW" "$CLEAN" "$IDXROOT" "$RUNBAM" "$QC" "$LOG" "$BASE/00_tables"

python3 - << 'PY' "$CSV" > "$BASE/00_tables/12JGI_run_rows.tsv"
import csv, sys
from pathlib import Path

csv_file = Path(sys.argv[1])
with open(csv_file, newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        species = row["Species"].strip()
        run = row["Run_Accession"].strip()
        layout = row["LibraryLayout"].strip()
        genome = row["GenomePath"].strip()
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

    python3 - << PY2 "$STDERR_FILE" "$RATE_FILE"
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
PY2

    echo "[$(date)] Done $SPECIES / $RUN"
done < "$BASE/00_tables/12JGI_run_rows.tsv"
