#!/bin/bash

# ======================================================================
# Script : 07_extract_rmod_log_summary.sh
# Purpose: Parse RepeatModeler rmod.log files across all genomes to
#          extract per-genome statistics reflecting inter-genome
#          variability (contig count, genome size, iterative round
#          count, repeat family discovery, classifier warnings, and
#          computational cost).
#
# Notes:
#   - The number of RepeatModeler sampling rounds is NOT fixed across
#     genomes; some genomes require additional rounds until no new
#     repeat families are identified. This script dynamically detects
#     the actual number of rounds run per genome rather than assuming
#     a fixed count.
#   - Program Time reported in rmod.log is wall-clock time for a
#     single job running with N threads (see THREADS_PER_JOB). CPU-core
#     hours are computed as wall-clock time x THREADS_PER_JOB to reflect
#     actual computational resource consumption.
# ======================================================================

set -e

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <repeat_working_dir> <output_csv_path> <threads_per_job>"
    echo "Example: $0 /path/to/04_repeat_working_dir /path/to/07_rmod_log_summary.csv 6"
    exit 1
fi

SRC_DIR=$(realpath "$1")
OUT_CSV="$2"
THREADS_PER_JOB="$3"

echo "Species,Sequences,Bases_bp,Rounds_Run,Families_Per_Round,Final_Round_Coverage_Percent,Total_Families_Discovered,Classifier_Warning,Program_Time_Wallclock,CPU_Core_Hours" > "$OUT_CSV"

for genome_dir in "$SRC_DIR"/*/; do
    genome_name=$(basename "$genome_dir")
    log_file=$(find "$genome_dir" -maxdepth 1 -name "*-rmod.log")

    if [ -z "$log_file" ]; then
        echo "$genome_name,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> "$OUT_CSV"
        echo "[MISSING] $genome_name has no -rmod.log file."
        continue
    fi

    sequences=$(grep -m1 "Sequences =" "$log_file" | sed -E 's/.*Sequences = ([0-9]+).*/\1/')
    bases=$(grep -m1 "Bases =" "$log_file" | sed -E 's/.*Bases = ([0-9]+).*/\1/')

    rounds_run=$(grep -c "RepeatModeler Round #" "$log_file")

    families_per_round=$(grep "families discovered" "$log_file" | sed -E 's/.*: ([0-9]+) families discovered.*/\1/' | paste -sd ";")

    final_coverage=$(grep "Input Database Coverage" "$log_file" | tail -1 | sed -E 's/.*\(\s*([0-9.]+)%\s*\).*/\1/')

    total_families=$(grep -m1 "discovery complete" "$log_file" | sed -E 's/.*: ([0-9]+) families found.*/\1/')

    if grep -q "WARNING" "$log_file"; then
        warning="Yes"
    else
        warning="No"
    fi

    runtime_str=$(grep -m1 "Program Time" "$log_file" | sed -E 's/.*Program Time: ([0-9]+):([0-9]+):([0-9]+).*/\1 \2 \3/')
    h=$(echo "$runtime_str" | awk '{print $1}')
    m=$(echo "$runtime_str" | awk '{print $2}')
    s=$(echo "$runtime_str" | awk '{print $3}')
    total_hours=$(echo "$h $m $s" | awk '{printf "%.4f", $1 + $2/60 + $3/3600}')
    cpu_hours=$(echo "$total_hours $THREADS_PER_JOB" | awk '{printf "%.4f", $1 * $2}')

    runtime_display="${h}:${m}:${s}"

    echo "$genome_name,$sequences,$bases,$rounds_run,\"$families_per_round\",$final_coverage,$total_families,$warning,$runtime_display,$cpu_hours" >> "$OUT_CSV"
done

echo "====================================================="
echo "Summary complete. Output saved to: $OUT_CSV"
echo "Total genomes processed: $(ls "$SRC_DIR" | wc -l)"
echo ""
echo "Distribution of Rounds_Run across genomes:"
tail -n +2 "$OUT_CSV" | awk -F',' '{print $4}' | sort | uniq -c
echo ""
echo "Total CPU-core-hours across full cohort:"
tail -n +2 "$OUT_CSV" | awk -F',' '{sum += $NF} END {printf "%.2f CPU-core-hours\n", sum}'
echo "====================================================="
