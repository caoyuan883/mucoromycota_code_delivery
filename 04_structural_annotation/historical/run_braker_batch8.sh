#!/usr/bin/env bash
set -euo pipefail

BASE="/media/desk16/tly8529/mucoromycota_annotation_braker"
GENOME_DIR="$BASE/01_Glomeromycotina_31"
PROT="$BASE/Fungi.fa"
OUTBASE="$BASE/04_braker_run"
LOGDIR="$OUTBASE/logs"

BRAKER="/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/BRAKER_src/scripts/braker.pl"

export GENEMARK_PATH="/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/aux/GeneMark/gmes_linux_64_4"
export PROTHINT_PATH="/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/aux/GeneMark/gmes_linux_64_4/ProtHint/bin"
export AUGUSTUS_CONFIG_PATH="/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/AUGUSTUS_src/config"
export AUGUSTUS_BIN_PATH="/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/AUGUSTUS_src/bin"
export AUGUSTUS_SCRIPTS_PATH="/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/AUGUSTUS_src/scripts"
export PYTHON3_PATH="$(which python3)"
export PERL5LIB="${PERL5LIB:-}"

mkdir -p "$OUTBASE" "$LOGDIR"
[ -f "$BRAKER" ] || { echo "ERROR: BRAKER not found: $BRAKER" >&2; exit 1; }
[ -f "$PROT" ] || { echo "ERROR: protein fasta not found: $PROT" >&2; exit 1; }
[ -d "$GENOME_DIR" ] || { echo "ERROR: genome dir not found: $GENOME_DIR" >&2; exit 1; }
[ -d "$GENEMARK_PATH" ] || { echo "ERROR: GENEMARK_PATH not found: $GENEMARK_PATH" >&2; exit 1; }
[ -d "$PROTHINT_PATH" ] || { echo "ERROR: PROTHINT_PATH not found: $PROTHINT_PATH" >&2; exit 1; }
[ -d "$AUGUSTUS_BIN_PATH" ] || { echo "ERROR: AUGUSTUS_BIN_PATH not found: $AUGUSTUS_BIN_PATH" >&2; exit 1; }
[ -d "$AUGUSTUS_SCRIPTS_PATH" ] || { echo "ERROR: AUGUSTUS_SCRIPTS_PATH not found: $AUGUSTUS_SCRIPTS_PATH" >&2; exit 1; }
[ -d "$AUGUSTUS_CONFIG_PATH" ] || { echo "ERROR: AUGUSTUS_CONFIG_PATH not found: $AUGUSTUS_CONFIG_PATH" >&2; exit 1; }

mapfile -t GENOMES < <(find "$GENOME_DIR" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | sort)

TOTAL=${#GENOMES[@]}
BATCH_SIZE=8
THREADS=8

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: no genome fasta files found in $GENOME_DIR" >&2
    exit 1
fi

echo "Start time: $(date)"
echo "Genome dir: $GENOME_DIR"
echo "Protein fasta: $PROT"
echo "Output dir: $OUTBASE"
echo "Total genomes: $TOTAL"
echo "Batch size: $BATCH_SIZE"
echo "Threads per sample: $THREADS"
echo
run_one () {
    local GENOME="$1"
    local NAME
    NAME=$(basename "$GENOME")
    NAME=${NAME%.fna}
    NAME=${NAME%.fa}
    NAME=${NAME%.fasta}

    local WD="$OUTBASE/$NAME"
    mkdir -p "$WD"

    echo "[START] $NAME  $(date)"
    nohup perl "$BRAKER" \
      --genome="$GENOME" \
      --prot_seq="$PROT" \
      --species="$NAME" \
      --workingdir="$WD" \
      --threads="$THREADS" \
      --gff3 \
      --softmasking \
      --skip_fixing_broken_genes \
      > "$LOGDIR/${NAME}.braker.log" 2>&1 &
}
i=0
batch=1

while [ $i -lt $TOTAL ]; do
    echo "========== BATCH $batch =========="
    for ((j=0; j<BATCH_SIZE && i<TOTAL; j++, i++)); do
        run_one "${GENOMES[$i]}"
        sleep 3
    done

    echo
    echo "[INFO] Waiting current batch to finish..."
    wait
    echo "[INFO] Batch $batch finished at $(date)"
    echo
    batch=$((batch+1))
done

echo "All batches finished at $(date)"

