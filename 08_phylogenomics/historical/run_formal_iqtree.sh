#!/usr/bin/env bash

ROOT="/media/desk16/tly8529/Mucoromycota_212_species_tree/IQTREE_run"

IQTREE="/media/desk16/tly8529/iqtree3_busco70/bin/iqtree-3.1.0-Linux/bin/iqtree3_intel"

MATRIX="$ROOT/04_supermatrix/Mucoromycota_212_300BUSCO_supermatrix.faa"

OUTDIR="$ROOT/06_formal_tree"

PREFIX="$OUTDIR/Mucoromycota_212_300BUSCO_JTTFG4"

mkdir -p "$OUTDIR"

echo "=================================================="
echo "FORMAL IQ-TREE ANALYSIS"
echo "=================================================="
echo "Date      : $(date)"
echo "IQ-TREE   : $IQTREE"
echo "Alignment : $MATRIX"
echo "Model     : JTT+F+G4"
echo "Threads   : 19"
echo "UFBoot    : 1000"
echo "SH-aLRT   : 1000"
echo "Seed      : 20260817"
echo "Prefix    : $PREFIX"
echo

if [[ ! -x "$IQTREE" ]]; then
    echo "ERROR: IQ-TREE executable was not found."
    return 1 2>/dev/null || true
fi

if [[ ! -s "$MATRIX" ]]; then
    echo "ERROR: input matrix was not found or is empty."
    return 1 2>/dev/null || true
fi

if [[ -e "${PREFIX}.treefile" ]] || [[ -e "${PREFIX}.iqtree" ]]; then
    echo "NOTICE: completed-result files already exist."
    echo "The formal analysis was not started again."
    return 0 2>/dev/null || true
fi

echo "Input check: PASS"
echo "Starting IQ-TREE..."

exec "$IQTREE" \
    -s "$MATRIX" \
    -st AA \
    -m JTT+F+G4 \
    -B 1000 \
    --alrt 1000 \
    -wbt \
    -alninfo \
    -T 19 \
    -mem 180G \
    -seed 20260817 \
    --prefix "$PREFIX"
