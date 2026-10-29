#!/bin/bash
# 由 03 调用：对一个 fasta 片段运行 InterProScan，带超时
f=$1; b=${2:-${f%.faa}}; A=${3:-$APPL}
n=$(grep -c '^>' $f)
T=$(( n*${TMO_PER:-30} )); [ $T -lt ${TMO_MIN:-1800} ] && T=${TMO_MIN:-1800}
[ -n "$TMO_FIX" ] && T=$TMO_FIX
\rm -rf $b.tmp; mkdir -p $b.tmp
start=$(date +%s)
timeout -k 120 $T $IPS -i $f -appl $A -goterms -dp -f TSV -cpu ${CCPU:-4} -b $b -T $b.tmp > $b.log 2>&1
rc=$?
\rm -rf $b.tmp
echo "$rc $(( $(date +%s)-start )) $T $n" > $b.rc
