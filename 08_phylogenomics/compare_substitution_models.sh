#!/usr/bin/env bash
# Step 03b: cheap restricted model comparison on the P matrix.
#
# Rationale: a full ModelFinder run over 10 models would cost 7-15 h on this
# matrix (409,430 distinct patterns, ~43 min per parameter optimisation).
# Instead we optimise parameters for a few candidate exchange matrices on the
# SAME fixed parsimony starting tree (-n 0, same -seed), all with +F+G4 so the
# number of free parameters is identical and log-likelihood ordering equals
# BIC ordering.
#
# LG+F+G4 was already computed by the memcheck stage (memcheck_P.log).
#
# Run inside tmux:
#   cd ~/dual_tree
#   bash compare_substitution_models.sh
#
# Takes roughly 45 min per model. Safe to re-run: finished models are skipped.

set -euo pipefail

# [release refactor] site paths come from the environment; all IQ-TREE options unchanged.
# historical: IQ=/media/desk16/tly8529/iqtree3_busco70/bin/iqtree-3.1.0-Linux/bin/iqtree3_intel
# historical: WORK=/media/desk16/tly8529/dual_tree
IQ=${IQ:?set IQ (path to the IQ-TREE 3 binary used for the logged runs)}
WORK=${WORK:?set WORK (folder with supermatrix_P_TAG.fasta and supermatrix_G_TAG.fasta)}
TAG=occ80_g90
THREADS=20
MEM=80G
SEED=12345
MP=$WORK/supermatrix_P_${TAG}.fasta

MODELS=("JTT+F+G4" "WAG+F+G4" "Q.pfam+F+G4")

cd "$WORK"
[ -x "$IQ" ] || { echo "ERROR: binary not executable: $IQ"; exit 1; }
[ -s "$MP" ] || { echo "ERROR: missing $MP"; exit 1; }

for M in "${MODELS[@]}"; do
  SAFE=$(echo "$M" | tr '+.' '__')
  PRE="mc_${SAFE}"
  if [ -s "${PRE}.log" ] && grep -q "BEST SCORE FOUND" "${PRE}.log"; then
    echo "== $M already done, skipping"
    continue
  fi
  echo "== $M  ($(date '+%F %T'))"
  "$IQ" -s "$MP" -st AA -m "$M" -n 0 \
        -T "$THREADS" -mem "$MEM" -seed "$SEED" --prefix "$PRE"
done

echo
echo "================ model comparison (P matrix, fixed parsimony tree) ================"
printf "%-16s %-22s %s\n" "MODEL" "LOG-LIKELIHOOD" "SOURCE"

LGLL=$(grep -h "BEST SCORE FOUND" memcheck_P.log 2>/dev/null | tail -1 | awk '{print $NF}')
[ -n "${LGLL:-}" ] && printf "%-16s %-22s %s\n" "LG+F+G4" "$LGLL" "memcheck_P.log"

for M in "${MODELS[@]}"; do
  SAFE=$(echo "$M" | tr '+.' '__')
  PRE="mc_${SAFE}"
  LL=$(grep -h "BEST SCORE FOUND" "${PRE}.log" 2>/dev/null | tail -1 | awk '{print $NF}')
  [ -n "${LL:-}" ] && printf "%-16s %-22s %s\n" "$M" "$LL" "${PRE}.log"
done

echo
echo "Pick the model with the HIGHEST (least negative) log-likelihood."
echo "All candidates share +F+G4, so the parameter count is identical."
echo "Then build both trees with that single model:"
echo "  bash run_iqtree_occ80.sh tree <MODEL>"
