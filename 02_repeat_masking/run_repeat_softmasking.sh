#!/bin/bash

# ======================================================================
# Script : run_repeat_softmasking.sh
# Purpose: Perform de novo repeat modeling and soft-masking on genome
#          assemblies using RepeatModeler and RepeatMasker in parallel.
#
# Dependencies (verified versions, extracted from run logs):
#   - RepeatModeler   v2.0.7
#   - Search Engine   rmblast 2.14.1+
#   - RECON           (bundled with RepeatModeler 2.0.7)
#   - RepeatScout     v1.0.7
#   - TRF             v4.10
#   - RepeatMasker    v4.2.3
#
# Note: LTR structural analysis module (-LTRStruct) is NOT enabled in
#       this pipeline. LTR retrotransposons are captured only through
#       the default de novo modules (RECON + RepeatScout), without
#       dedicated structural validation.
# ======================================================================

set -e

# Help message
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_fasta_dir> <output_dir> <threads_per_job>"
    echo "Example: $0 /path/to/genomes /path/to/output_masking 6"
    exit 1
fi

INPUT_DIR=$(realpath "$1")
BASE_OUT=$(realpath "$2")
THREADS_PER_JOB="$3"
# Number of parallel genomes to process simultaneously
PARALLEL_JOBS=10

# Setup output directories
OUT_MASKED="$BASE_OUT/01_softmasked_genomes"
OUT_RM_WORK="$BASE_OUT/02_repeat_working_dir"
LOG_DIR="$BASE_OUT/03_parallel_logs"

mkdir -p "$OUT_MASKED" "$OUT_RM_WORK" "$LOG_DIR"

# Check dependencies
for cmd in BuildDatabase RepeatModeler RepeatMasker; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd could not be found. Please ensure it is in your PATH."
        exit 1
    fi
done

# Define the processing function
process_genome() {
    GENOME="$1"
    FILENAME=$(basename "$GENOME")
    PREFIX="${FILENAME%.*}"
    LOG_FILE="$LOG_DIR/${PREFIX}.log"

    if [ -f "$OUT_MASKED/${PREFIX}_softmasked.fa" ]; then
        echo "[SKIP] Soft-masked file for $PREFIX already exists."
        return 0
    fi

    echo "=====================================================" > "$LOG_FILE"
    echo "Processing Genome: $PREFIX" >> "$LOG_FILE"
    date >> "$LOG_FILE"

    WORK_DIR="$OUT_RM_WORK/$PREFIX"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR" || return 1

    echo "[1/3] Running BuildDatabase..." >> "$LOG_FILE"
    BuildDatabase -name "$PREFIX" "$GENOME" >> "$LOG_FILE" 2>&1

    echo "[2/3] Running RepeatModeler (v2.0.7; engine: rmblast 2.14.1+; LTRStruct disabled)..." >> "$LOG_FILE"
    RepeatModeler -database "$PREFIX" -threads "$THREADS_PER_JOB" >> "$LOG_FILE" 2>&1

    if [ -f "${PREFIX}-families.fa" ]; then
        echo "[3/3] Running RepeatMasker (v4.2.3; -xsmall soft-masking)..." >> "$LOG_FILE"
        RepeatMasker -lib "${PREFIX}-families.fa" -pa "$THREADS_PER_JOB" -xsmall -gff -dir "$WORK_DIR" "$GENOME" >> "$LOG_FILE" 2>&1

        if [ -f "${FILENAME}.masked" ]; then
            cp "${FILENAME}.masked" "$OUT_MASKED/${PREFIX}_softmasked.fa"
            echo "[SUCCESS] Masking completed for $PREFIX" >> "$LOG_FILE"
        else
            echo "[ERROR] .masked file not generated for $PREFIX. Check logs." >> "$LOG_FILE"
        fi
    else
        echo "[ERROR] RepeatModeler family library not generated for $PREFIX. Check logs." >> "$LOG_FILE"
    fi
    date >> "$LOG_FILE"
}

export -f process_genome
export OUT_MASKED OUT_RM_WORK LOG_DIR THREADS_PER_JOB

echo "Starting parallel soft-masking for genomes in: $INPUT_DIR"
echo "Results will be saved to: $BASE_OUT"

# Execute in parallel
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.gna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | \
    xargs -I {} -P "$PARALLEL_JOBS" bash -c 'process_genome "$@"' _ {}

echo "All genome masking tasks have been submitted!"
