#!/bin/bash
# ============================================================================
# Step 09 — 重建权威 manifest（只读，不动任何产物）
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

FULL_ROOT=$BATCH_ROOT/_bbc_full
OUT=$BATCH_ROOT/_manifest_v4.tsv
CONS=$BATCH_ROOT/_busco_consistency.tsv
REP=$BATCH_ROOT/_manifest_v4_report.txt

getv() { local k=$1 f=$2; grep "^${k}:" "$f" 2>/dev/null | head -1 | sed "s/^${k}://; s/%.*//; s/^ *//"; }
nz()   { [[ -z "$1" ]] && echo "NA" || echo "$1"; }

echo -e "species\tn_genes\tn_proteins\tbusco_C\tbusco_S\tbusco_D\tbusco_F\tbusco_M\tbusco_n\tnative_braker_M\tnative_augustus_M\tnative_genemark_M\tnative_decision\tverify_braker_M\tverify_augustus_M\tverify_genemark_M\tverify_mode\tfinal_gene_set\tconsistency" > $OUT
echo -e "species\tset\tnative_M\tverify_M\tdelta\tflag" > $CONS

n=0
while read sp; do
    [[ -z "$sp" ]] && continue
    n=$((n+1))
    wd=$BATCH_ROOT/$sp
    vf=$FULL_ROOT/$sp

    n_prot=$(grep -c '^>' $wd/braker.aa 2>/dev/null); n_prot=$(nz "$n_prot")
    n_gene=$(awk -F'\t' '$3=="gene"' $wd/braker.gtf 2>/dev/null | wc -l)
    if [[ "$n_gene" == "0" ]]; then
        n_gene=$(awk -F'\t' '$3=="transcript" || $3=="mRNA"' $wd/braker.gtf 2>/dev/null | wc -l)
    fi
    [[ "$n_gene" == "0" ]] && n_gene="NA"

    nb=$(nz "$(getv M $wd/bbc/braker/summary.txt)")
    na=$(nz "$(getv M $wd/bbc/augustus/summary.txt)")
    ng=$(nz "$(getv M $wd/bbc/genemark/summary.txt)")

    nd="NA"
    if [[ -s "$wd/best_by_compleasm.log" ]]; then
        if grep -q "is the best one" $wd/best_by_compleasm.log; then
            nd="keep_braker_best"
        elif grep -q "there are no BUSCOs to be added" $wd/best_by_compleasm.log; then
            nd="keep_braker_no_add"
        elif grep -qiE "will be added|merge" $wd/best_by_compleasm.log; then
            nd="merged"
        else
            nd="unparsed"
        fi
    fi

    vb=$(nz "$(getv M $vf/braker/summary.txt)")
    va=$(nz "$(getv M $vf/augustus/summary.txt)")
    vg=$(nz "$(getv M $vf/genemark/summary.txt)")
    vC=$(nz "$(getv C $vf/braker/summary.txt)")
    vS=$(nz "$(getv S $vf/braker/summary.txt)")
    vD=$(nz "$(getv D $vf/braker/summary.txt)")
    vF=$(nz "$(getv F $vf/braker/summary.txt)")
    vN=$(grep '^N:' $vf/braker/summary.txt 2>/dev/null | head -1 | sed 's/^N://; s/^ *//')
    vN=$(nz "$vN")
    if [[ "$vC" == "NA" && "$vS" != "NA" && "$vD" != "NA" ]]; then
        vC=$(awk -v a=$vS -v b=$vD 'BEGIN{printf "%.2f", a+b}')
    fi

    mode="full"
    [[ "$va" == "NA" ]] && mode="no_augustus"

    fset="braker.gtf"
    [[ -s "$wd/braker_final.gtf" ]] && fset="braker_final.gtf"

    cons="OK"
    for pair in "braker:$nb:$vb" "augustus:$na:$va" "genemark:$ng:$vg"; do
        setn=${pair%%:*}; rest=${pair#*:}; nv=${rest%%:*}; vv=${rest##*:}
        if [[ "$nv" == "NA" || "$vv" == "NA" ]]; then
            echo -e "$sp\t$setn\t$nv\t$vv\tNA\tNOT_COMPARABLE" >> $CONS
            continue
        fi
        d=$(awk -v a=$nv -v b=$vv 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.2f", d}')
        fl="OK"
        awk -v d=$d 'BEGIN{exit !(d>0.5)}' && { fl="DIVERGENT"; cons="DIVERGENT"; }
        echo -e "$sp\t$setn\t$nv\t$vv\t$d\t$fl" >> $CONS
    done

    echo -e "$sp\t$n_gene\t$n_prot\t$vC\t$vS\t$vD\t$vF\t$vb\t$vN\t$nb\t$na\t$ng\t$nd\t$vb\t$va\t$vg\t$mode\t$fset\t$cons" >> $OUT
    printf "\r  处理中 %3d/127  %-40s" "$n" "$sp"
done < $SPECIES_LIST_FILE
echo ""

{
echo "=========================================================="
echo " BRAKER2 v3 — Mucoromycota 注释资源 权威汇总 (manifest v4)"
echo " 生成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================================="
echo ""
echo "物种总数: $(awk 'NR>1' $OUT | wc -l)"
echo ""
echo "--- 最终基因集 ---"
awk -F'\t' 'NR>1 {c[$18]++} END {for(k in c) printf "  %-20s %d\n", k, c[k]}' $OUT
echo ""
echo "--- BRAKER 原生决策 ---"
awk -F'\t' 'NR>1 {c[$13]++} END {for(k in c) printf "  %-22s %d\n", k, c[k]}' $OUT
echo ""
echo "--- 复现模式 ---"
awk -F'\t' 'NR>1 {c[$17]++} END {for(k in c) printf "  %-16s %d\n", k, c[k]}' $OUT
echo ""
echo "--- BUSCO 完整度 C% (稳定库复现) ---"
awk -F'\t' 'NR>1 && $4!="NA" {s+=$4; n++; if(min==""||$4<min)min=$4; if($4>max)max=$4} END {printf "  n=%d  mean=%.2f  min=%.2f  max=%.2f\n", n, s/n, min, max}' $OUT
echo ""
echo "  C >= 90%:  $(awk -F'\t' 'NR>1 && $4!="NA" && $4>=90' $OUT | wc -l)"
echo "  C 80-90%:  $(awk -F'\t' 'NR>1 && $4!="NA" && $4>=80 && $4<90' $OUT | wc -l)"
echo "  C 70-80%:  $(awk -F'\t' 'NR>1 && $4!="NA" && $4>=70 && $4<80' $OUT | wc -l)"
echo "  C <  70%:  $(awk -F'\t' 'NR>1 && $4!="NA" && $4<70' $OUT | wc -l)"
echo ""
echo "--- 蛋白数 ---"
awk -F'\t' 'NR>1 && $3!="NA" {s+=$3; n++; if(min==""||$3<min)min=$3; if($3>max)max=$3} END {printf "  n=%d  mean=%.0f  min=%d  max=%d  total=%d\n", n, s/n, min, max, s}' $OUT
echo ""
echo "--- 两轮打分一致性 ---"
echo "  可比对条目: $(awk -F'\t' 'NR>1 && $6!="NOT_COMPARABLE"' $CONS | wc -l)"
echo "  一致 (|Δ|<=0.5): $(awk -F'\t' 'NR>1 && $6=="OK"' $CONS | wc -l)"
echo "  分歧 (|Δ|> 0.5): $(awk -F'\t' 'NR>1 && $6=="DIVERGENT"' $CONS | wc -l)"
echo "  不可比 (缺一侧): $(awk -F'\t' 'NR>1 && $6=="NOT_COMPARABLE"' $CONS | wc -l)"
echo ""
echo "  分歧最大的 10 条:"
awk -F'\t' 'NR>1 && $6=="DIVERGENT"' $CONS | sort -t$'\t' -k5,5gr | head -10 | awk -F'\t' '{printf "    %-34s %-9s native=%-6s verify=%-6s Δ=%s\n", $1, $2, $3, $4, $5}'
echo ""
echo "  受影响物种数: $(awk -F'\t' 'NR>1 && $19=="DIVERGENT"' $OUT | wc -l) / 127"
echo ""
echo "--- 缺失最严重的 10 株 (verify braker M%) ---"
awk -F'\t' 'NR>1 && $8!="NA" {printf "  %-36s M=%-6s C=%s\n", $1, $8, $4}' $OUT | sort -t= -k2,2gr | head -10
echo ""
echo "=========================================================="
} > $REP

cat $REP
echo ""
echo "输出文件:"
echo "  $OUT"
echo "  $CONS"
echo "  $REP"
