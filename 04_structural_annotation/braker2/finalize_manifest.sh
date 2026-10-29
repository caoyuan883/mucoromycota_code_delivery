#!/bin/bash
# ============================================================================
# Step 09b — 最终权威 manifest（合并 BRAKER 原生决策 + 06c 补跑决策）
# 只读，不修改任何注释产物
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

FULL_ROOT=$BATCH_ROOT/_bbc_full
OUT=$BATCH_ROOT/_manifest_final.tsv
REP=$BATCH_ROOT/_manifest_final_report.txt

getv() { grep "^${1}:" "$2" 2>/dev/null | head -1 | sed "s/^${1}://; s/%.*//; s/^ *//"; }
nz()   { [[ -z "$1" ]] && echo "NA" || echo "$1"; }

parse_dec() {
    local f=$1
    [[ -s "$f" ]] || { echo "NA"; return; }
    if   grep -q "is the best one" "$f";                then echo "keep_braker_best"
    elif grep -q "there are no BUSCOs to be added" "$f"; then echo "keep_braker_no_add"
    elif grep -qiE "will be added to the braker|merged" "$f"; then echo "merged"
    elif grep -q "Error in file" "$f";                   then echo "ERROR"
    else echo "unparsed"; fi
}

echo -e "species\tn_genes\tn_proteins\tbusco_C\tbusco_S\tbusco_D\tbusco_F\tbusco_M\tbusco_n\tbraker_M\taugustus_M\tgenemark_M\tdecision\tdecision_source\tbraker_native_bbc\tfinal_gene_set" > $OUT

n=0
while read sp; do
    [[ -z "$sp" ]] && continue
    n=$((n+1))
    wd=$BATCH_ROOT/$sp
    vf=$FULL_ROOT/$sp

    n_prot=$(nz "$(grep -c '^>' $wd/braker.aa 2>/dev/null)")
    n_gene=$(awk -F'\t' '$3=="gene"' $wd/braker.gtf 2>/dev/null | wc -l)
    [[ "$n_gene" == "0" ]] && n_gene=$(awk -F'\t' '$3=="transcript" || $3=="mRNA"' $wd/braker.gtf 2>/dev/null | wc -l)
    [[ "$n_gene" == "0" ]] && n_gene="NA"

    nat_dec=$(parse_dec "$wd/best_by_compleasm.log")
    pos_dec=$(parse_dec "$vf.log")

    if [[ "$nat_dec" == "keep_braker_best" || "$nat_dec" == "keep_braker_no_add" || "$nat_dec" == "merged" ]]; then
        dec=$nat_dec; src="braker_native"
        bM=$(nz "$(getv M $wd/bbc/braker/summary.txt)")
        aM=$(nz "$(getv M $wd/bbc/augustus/summary.txt)")
        gM=$(nz "$(getv M $wd/bbc/genemark/summary.txt)")
        C=$(nz "$(getv C $wd/bbc/braker/summary.txt)")
        S=$(nz "$(getv S $wd/bbc/braker/summary.txt)")
        D=$(nz "$(getv D $wd/bbc/braker/summary.txt)")
        F=$(nz "$(getv F $wd/bbc/braker/summary.txt)")
        N=$(nz "$(grep '^N:' $wd/bbc/braker/summary.txt 2>/dev/null | head -1 | sed 's/^N://; s/^ *//')")
        nat="yes"
    else
        dec=$pos_dec; src="post_hoc_06c"
        bM=$(nz "$(getv M $vf/braker/summary.txt)")
        aM=$(nz "$(getv M $vf/augustus/summary.txt)")
        gM=$(nz "$(getv M $vf/genemark/summary.txt)")
        C=$(nz "$(getv C $vf/braker/summary.txt)")
        S=$(nz "$(getv S $vf/braker/summary.txt)")
        D=$(nz "$(getv D $vf/braker/summary.txt)")
        F=$(nz "$(getv F $vf/braker/summary.txt)")
        N=$(nz "$(grep '^N:' $vf/braker/summary.txt 2>/dev/null | head -1 | sed 's/^N://; s/^ *//')")
        nat="no_braker_internal_failed"
    fi

    if [[ "$C" == "NA" && "$S" != "NA" && "$D" != "NA" ]]; then
        C=$(awk -v a=$S -v b=$D 'BEGIN{printf "%.2f", a+b}')
    fi

    fset="braker.gtf"
    [[ -s "$wd/braker_final.gtf" ]] && fset="braker_final.gtf"

    echo -e "$sp\t$n_gene\t$n_prot\t$C\t$S\t$D\t$F\t$bM\t$N\t$bM\t$aM\t$gM\t$dec\t$src\t$nat\t$fset" >> $OUT
    printf "\r  %3d/127  %-42s" "$n" "$sp"
done < $SPECIES_LIST_FILE
echo ""

{
echo "=========================================================="
echo " Mucoromycota BRAKER2 注释资源 — 最终 manifest"
echo " 生成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================================="
echo ""
echo "物种总数: $(awk 'NR>1' $OUT | wc -l)"
echo ""
echo "--- 择优决策（全部 127 株）---"
awk -F'\t' 'NR>1 {c[$13]++} END {for(k in c) printf "  %-24s %d\n", k, c[k]}' $OUT
echo ""
echo "--- 决策来源 ---"
awk -F'\t' 'NR>1 {c[$14]++} END {for(k in c) printf "  %-18s %d\n", k, c[k]}' $OUT
echo ""
echo "--- 最终基因集 ---"
awk -F'\t' 'NR>1 {c[$16]++} END {for(k in c) printf "  %-20s %d\n", k, c[k]}' $OUT
echo ""
echo "--- BUSCO 完整度 C% ---"
awk -F'\t' 'NR>1 && $4!="NA" {s+=$4; n++; if(min==""||$4<min)min=$4; if($4>max)max=$4} END {printf "  n=%d  mean=%.2f  min=%.2f  max=%.2f\n", n, s/n, min, max}' $OUT
echo "  C>=95%: $(awk -F'\t' 'NR>1 && $4!="NA" && $4>=95' $OUT | wc -l)   C>=90%: $(awk -F'\t' 'NR>1 && $4!="NA" && $4>=90' $OUT | wc -l)   C<80%: $(awk -F'\t' 'NR>1 && $4!="NA" && $4<80' $OUT | wc -l)"
echo ""
echo "--- 蛋白数 ---"
awk -F'\t' 'NR>1 && $3!="NA" {s+=$3; n++; if(min==""||$3<min)min=$3; if($3>max)max=$3} END {printf "  n=%d  mean=%.0f  min=%d  max=%d  total=%d\n", n, s/n, min, max, s}' $OUT
echo ""
echo "--- 异常检查 ---"
echo "  decision=ERROR:    $(awk -F'\t' 'NR>1 && $13=="ERROR"' $OUT | wc -l)"
echo "  decision=unparsed: $(awk -F'\t' 'NR>1 && $13=="unparsed"' $OUT | wc -l)"
echo "  decision=NA:       $(awk -F'\t' 'NR>1 && $13=="NA"' $OUT | wc -l)"
echo "  BUSCO 缺值:        $(awk -F'\t' 'NR>1 && $4=="NA"' $OUT | wc -l)"
echo "  基因数缺值:        $(awk -F'\t' 'NR>1 && $2=="NA"' $OUT | wc -l)"
awk -F'\t' 'NR>1 && ($13=="ERROR"||$13=="unparsed"||$13=="NA"||$4=="NA"||$2=="NA") {print "    异常: "$1"  dec="$13"  C="$4"  genes="$2}' $OUT
echo ""
echo "=========================================================="
} > $REP

cat $REP
echo ""
echo "输出: $OUT"
