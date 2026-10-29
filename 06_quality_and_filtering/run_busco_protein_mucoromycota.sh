#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script : run_busco_protein_mucoromycota.sh
# Purpose: BUSCO protein-mode assessment (mucoromycota_odb12 only, offline)
#          for *_protein.faa files.
#
# Example:
#   bash run_busco_protein_mucoromycota.sh \
#     --protein-dir /home/py/Desktop/147_protein \
#     --out-dir /home/py/Desktop/147_protein/busco_protein \
#     --db-path /path/to/busco_db \
#     --threads 16
# =============================================================================

PROTEIN_DIR=""
OUT_DIR=""
DB_PATH=""
THREADS=16
CONDA_ENV="busco6"
LINEAGE="mucoromycota_odb12"
FORCE_RERUN=0

usage() {
cat <<'USAGE_EOF'
Usage:
  bash run_busco_protein_mucoromycota.sh \
    --protein-dir /path/to/proteins \
    --out-dir /path/to/busco_output \
    --db-path /path/to/busco_db \
    [--threads 16] [--conda-env busco6] [--lineage mucoromycota_odb12] [--force-rerun]
USAGE_EOF
exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --protein-dir) PROTEIN_DIR="$2"; shift 2 ;;
    --out-dir)     OUT_DIR="$2"; shift 2 ;;
    --db-path)     DB_PATH="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --conda-env)   CONDA_ENV="$2"; shift 2 ;;
    --lineage)     LINEAGE="$2"; shift 2 ;;
    --force-rerun) FORCE_RERUN=1; shift 1 ;;
    --help|-h)     usage ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage ;;
  esac
done

if [[ -z "$PROTEIN_DIR" || -z "$OUT_DIR" || -z "$DB_PATH" ]]; then
  echo "ERROR: --protein-dir, --out-dir, and --db-path are required." >&2
  usage
fi
[[ -d "$PROTEIN_DIR" ]] || { echo "ERROR: Protein directory not found: $PROTEIN_DIR" >&2; exit 1; }
[[ -d "$DB_PATH" ]]     || { echo "ERROR: BUSCO database directory not found: $DB_PATH" >&2; exit 1; }

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

LOG_DIR="$OUT_DIR/logs"
RES_DIR="$OUT_DIR/results_mucoromycota"
MASTER_LOG="$LOG_DIR/busco_master.log"
mkdir -p "$LOG_DIR" "$RES_DIR"

cleanup() {
  echo "[$(date)] INTERRUPTED: stopping BUSCO batch" | tee -a "$MASTER_LOG"
  exit 130
}
trap cleanup INT TERM

echo "===== START $LINEAGE =====" | tee -a "$MASTER_LOG"

shopt -s nullglob
files=( "$PROTEIN_DIR"/*_protein.faa )
shopt -u nullglob
if [[ "${#files[@]}" -eq 0 ]]; then
  echo "[$(date)] ERROR: no *_protein.faa files found in $PROTEIN_DIR" | tee -a "$MASTER_LOG" >&2
  exit 1
fi
echo "[$(date)] Found ${#files[@]} proteome file(s)" | tee -a "$MASTER_LOG"

for f in "${files[@]}"; do
  base=$(basename "$f")
  sample="${base%_protein.faa}"

  sample_dir="$RES_DIR/$sample"
  sample_log="$LOG_DIR/${sample}.${LINEAGE}.log"

  if [[ -d "$sample_dir" ]] && ls "$sample_dir"/short_summary*.txt >/dev/null 2>&1; then
    echo "[$(date)] SKIP $sample (complete result exists)" | tee -a "$MASTER_LOG"
    continue
  fi

  if [[ -d "$sample_dir" ]]; then
    if [[ "$FORCE_RERUN" -eq 0 ]]; then
      echo "[$(date)] ERROR $sample (incomplete result directory exists; use --force-rerun)" | tee -a "$MASTER_LOG" >&2
      continue
    else
      rm -rf "$sample_dir"
    fi
  fi

  echo "[$(date)] START $sample" | tee -a "$MASTER_LOG"

  set +e
  busco \
    -i "$f" \
    -m protein \
    -l "$LINEAGE" \
    --offline \
    --download_path "$DB_PATH" \
    -o "$sample" \
    --out_path "$RES_DIR" \
    -c "$THREADS" \
    > "$sample_log" 2>&1
  status=$?
  set -e

  if [[ "$status" -ne 0 ]]; then
    echo "[$(date)] FAIL $sample (exit=$status)" | tee -a "$MASTER_LOG" >&2
  else
    echo "[$(date)] END $sample" | tee -a "$MASTER_LOG"
  fi
done

echo "===== END $LINEAGE =====" | tee -a "$MASTER_LOG"
echo "[$(date)] ALL BUSCO RUNS FINISHED" | tee -a "$MASTER_LOG"
