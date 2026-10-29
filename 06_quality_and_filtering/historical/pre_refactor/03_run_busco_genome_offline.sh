#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script : 03_run_busco_genome_offline.sh
# Purpose: Run BUSCO genome-mode assessment on genome FASTA files using
#          offline lineage datasets.
#
# Example:
#   bash 03_run_busco_genome_offline.sh \
#     --genome-dir /path/to/genomes \
#     --out-dir /path/to/busco_output \
#     --db-path /path/to/busco_db \
#     --threads 16
#
# Notes:
#   - Both mucoromycota_odb12 and fungi_odb12 are run for each genome.
#   - Complete existing results are skipped automatically.
#   - Incomplete existing result directories are NOT overwritten unless
#     --force-rerun is provided.
# =============================================================================

GENOME_DIR=""
OUT_DIR=""
DB_PATH=""
THREADS=16
CONDA_ENV="busco6"
FORCE_RERUN=0

usage() {
  cat <<'USAGE_EOF'
Usage:
  bash 03_run_busco_genome_offline.sh \
    --genome-dir /path/to/genomes \
    --out-dir /path/to/busco_output \
    --db-path /path/to/busco_db \
    [--threads 16] [--conda-env busco6] [--force-rerun]

Required:
  --genome-dir   Directory containing input genome FASTA files (*.fna)
  --out-dir      Output root directory for BUSCO results
  --db-path      BUSCO offline database directory

Optional:
  --threads      Number of CPU threads per BUSCO run (default: 16)
  --conda-env    Conda environment name (default: busco6)
  --force-rerun  Remove incomplete existing sample directories before rerun
  --help         Show this help message
USAGE_EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome-dir)
      GENOME_DIR="$2"; shift 2 ;;
    --out-dir)
      OUT_DIR="$2"; shift 2 ;;
    --db-path)
      DB_PATH="$2"; shift 2 ;;
    --threads)
      THREADS="$2"; shift 2 ;;
    --conda-env)
      CONDA_ENV="$2"; shift 2 ;;
    --force-rerun)
      FORCE_RERUN=1; shift 1 ;;
    --help|-h)
      usage ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage ;;
  esac
done

if [[ -z "$GENOME_DIR" || -z "$OUT_DIR" || -z "$DB_PATH" ]]; then
  echo "ERROR: --genome-dir, --out-dir, and --db-path are required." >&2
  usage
fi

if [[ ! -d "$GENOME_DIR" ]]; then
  echo "ERROR: Genome directory not found: $GENOME_DIR" >&2
  exit 1
fi

if [[ ! -d "$DB_PATH" ]]; then
  echo "ERROR: BUSCO database directory not found: $DB_PATH" >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

LOG_DIR="$OUT_DIR/logs"
OUT1="$OUT_DIR/results_mucoromycota"
OUT2="$OUT_DIR/results_fungi"
MASTER_LOG="$LOG_DIR/busco_master.log"

mkdir -p "$LOG_DIR" "$OUT1" "$OUT2"

cleanup() {
  echo "[$(date)] INTERRUPTED: stopping BUSCO batch" | tee -a "$MASTER_LOG"
  exit 130
}
trap cleanup INT TERM

run_lineage() {
  local lineage="$1"
  local outdir="$2"

  echo "===== START $lineage =====" | tee -a "$MASTER_LOG"

  for f in "$GENOME_DIR"/*.fna; do
    [[ -e "$f" ]] || continue

    base=$(basename "$f")
    sample="${base%.fna}"
    sample="${sample%_genomic}"

    sample_dir="$outdir/$sample"
    sample_log="$LOG_DIR/${sample}.${lineage}.log"

    if [[ -d "$sample_dir" ]] && ls "$sample_dir"/short_summary*.txt >/dev/null 2>&1; then
      echo "[$(date)] SKIP $sample $lineage (complete result exists)" | tee -a "$MASTER_LOG"
      continue
    fi

    if [[ -d "$sample_dir" ]]; then
      if [[ "$FORCE_RERUN" -eq 0 ]]; then
        echo "[$(date)] ERROR $sample $lineage (incomplete result directory exists; use --force-rerun to overwrite)" | tee -a "$MASTER_LOG" >&2
        return 1
      else
        rm -rf "$sample_dir"
      fi
    fi

    echo "[$(date)] START $sample $lineage" | tee -a "$MASTER_LOG"

    busco \
      -i "$f" \
      -m genome \
      -l "$lineage" \
      --offline \
      --download_path "$DB_PATH" \
      -o "$sample" \
      --out_path "$outdir" \
      -c "$THREADS" \
      > "$sample_log" 2>&1

    status=$?
    if [[ "$status" -ne 0 ]]; then
      echo "[$(date)] FAIL $sample $lineage (exit=$status)" | tee -a "$MASTER_LOG" >&2
    else
      echo "[$(date)] END $sample $lineage" | tee -a "$MASTER_LOG"
    fi
  done

  echo "===== END $lineage =====" | tee -a "$MASTER_LOG"
}

run_lineage "mucoromycota_odb12" "$OUT1"
run_lineage "fungi_odb12" "$OUT2"

echo "[$(date)] ALL BUSCO RUNS FINISHED" | tee -a "$MASTER_LOG"
