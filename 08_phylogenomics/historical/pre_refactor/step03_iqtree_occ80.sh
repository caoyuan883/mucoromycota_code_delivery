#!/usr/bin/env bash
# Step 03: IQ-TREE runs for the dual-tree validation (occ80, g90 matrices).
#
# Run this INSIDE tmux, one stage at a time:
#   tmux new -s dualtree
#   cd ~/dual_tree
#   bash step03_iqtree_occ80.sh memcheck
#   bash step03_iqtree_occ80.sh mf
#   bash step03_iqtree_occ80.sh tree JTT+F+G4        # model from the mf stage
#
# Same binary, same threads, same seed, same model for every run.
# Re-running any stage with the same command resumes from the IQ-TREE
# checkpoint; never add -redo.

set -euo pipefail

IQ=/media/desk16/tly8529/iqtree3_busco70/bin/iqtree-3.1.0-Linux/bin/iqtree3_intel
WORK=/media/desk16/tly8529/dual_tree
TAG=occ80_g90
THREADS=20
MEM=80G
SEED=12345
MP=$WORK/supermatrix_P_${TAG}.fasta
MG=$WORK/supermatrix_G_${TAG}.fasta

STAGE=${1:-}
MODEL=${2:-}

cd "$WORK"

check_env() {
  [ -x "$IQ" ] || { echo "ERROR: IQ-TREE binary not executable: $IQ"; exit 1; }
  [ -s "$MP" ] || { echo "ERROR: missing $MP"; exit 1; }
  [ -s "$MG" ] || { echo "ERROR: missing $MG"; exit 1; }
  echo "--- binary"
  "$IQ" --version | head -2
  echo "--- input"
  for f in "$MP" "$MG"; do
    printf "%s  seqs=%s  size=%s\n" "$(basename "$f")" \
      "$(grep -c '>' "$f")" "$(du -h "$f" | cut -f1)"
  done
  echo "--- machine"
  uptime
  free -g | head -2
  echo
}

case "$STAGE" in

  memcheck)
    check_env
    echo "=== stage: memcheck (parameters on a parsimony tree, no tree search) ==="
    "$IQ" -s "$MP" -st AA -m LG+F+G4 -n 0 \
          -T "$THREADS" -mem "$MEM" -seed "$SEED" --prefix memcheck_P
    echo
    echo "=== key lines ==="
    grep -iE "RAM|memory|distinct patterns|Threads|sequences with" memcheck_P.log || true
    echo
    echo "Next: if required RAM is comfortably below the free memory shown above,"
    echo "      run:  bash $0 mf"
    ;;

  mf)
    check_env
    echo "=== stage: restricted model selection on the P matrix (BIC) ==="
    "$IQ" -s "$MP" -st AA -m MF \
          -mset LG,JTT,WAG,Q.pfam,Q.yeast -mfreq F -mrate G,R \
          -cmin 4 -cmax 4 -merit BIC \
          -T "$THREADS" -mem "$MEM" -seed "$SEED" --prefix MF_P_${TAG}
    echo
    echo "=== key lines ==="
    grep -E "Best-fit model|BIC score" MF_P_${TAG}.iqtree | head -5 || true
    grep -E "Best-fit model" MF_P_${TAG}.log | head -3 || true
    echo
    echo "Next: take the model printed above and run:"
    echo "      bash $0 tree <MODEL>"
    ;;

  tree)
    [ -n "$MODEL" ] || { echo "ERROR: give the model, e.g. bash $0 tree JTT+F+G4"; exit 1; }
    check_env
    echo "=== stage: ML trees, model fixed to $MODEL for BOTH sides ==="
    for SIDE in P G; do
      MAT=$WORK/supermatrix_${SIDE}_${TAG}.fasta
      PRE=${SIDE}_${TAG}
      echo
      echo "--- $SIDE side -> $PRE"
      date
      "$IQ" -s "$MAT" -st AA -m "$MODEL" \
            -B 1000 -alrt 1000 \
            -T "$THREADS" -mem "$MEM" -seed "$SEED" --prefix "$PRE"
      date
    done
    echo
    echo "=== summary ==="
    for SIDE in P G; do
      PRE=${SIDE}_${TAG}
      echo "--- $PRE"
      grep -E "Model of substitution|Log-likelihood of the tree|Total wall-clock time" \
           ${PRE}.iqtree ${PRE}.log 2>/dev/null | head -4 || true
      [ -s ${PRE}.treefile ] && echo "treefile OK: $(wc -c < ${PRE}.treefile) bytes"
    done
    echo
    echo "Next: quick unrooted RF check including the 5 outgroups (rough only):"
    echo "  $IQ -rf P_${TAG}.treefile G_${TAG}.treefile"
    echo "The reported result must come from the 202-leaf pruned comparison (step 04)."
    ;;

  *)
    echo "usage: bash $0 {memcheck|mf|tree <MODEL>}"
    exit 1
    ;;
esac
