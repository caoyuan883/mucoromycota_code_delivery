#!/bin/bash
# ============================================================================
# Step 10 — 为 BRAKER 中断的 101 株补生成 braker.gff3
#  工具与参数和 BRAKER 原生调用一致（已用原生 24 株验证 diff=0）
#  安全: 已存在 gff3 的物种一律跳过，绝不覆盖
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

# [release refactor] historical: G2G=/media/desk16/tly8529/glo_braker_test/08_env/braker2_bundle_transfer/src/AUGUSTUS_src/scripts/gtf2gff.pl
G2G=$ENV_ROOT/braker2_bundle_transfer/src/AUGUSTUS_src/scripts/gtf2gff.pl
LOG=$BATCH_ROOT/_gff3_backfill.log

log_g() { echo "[$(ts)] [GFF3] $*" | tee -a $LOG; }

log_g "===== gff3 补生成 ====="
log_g "  工具: $G2G"
log_g "  参数: --gff3 --printExon --printIntron"

n=0; made=0; skip=0; fail=0
while read sp; do
    [[ -z "$sp" ]] && continue
    n=$((n+1))
    wd=$BATCH_ROOT/$sp

    if [[ -s "$wd/braker.gff3" ]]; then
        skip=$((skip+1)); continue
    fi
    if [[ ! -s "$wd/braker.gtf" ]]; then
        log_g "  [$sp] FAIL 无 braker.gtf"; fail=$((fail+1)); continue
    fi

    tmp=$wd/.braker.gff3.tmp
    cat $wd/braker.gtf | perl $G2G --gff3 --printExon --printIntron --out=$tmp 2> $wd/errors/gtf2gff.backfill.stderr
    rc=$?

    if [[ $rc -ne 0 || ! -s "$tmp" ]]; then
        log_g "  [$sp] FAIL rc=$rc  $(head -1 $wd/errors/gtf2gff.backfill.stderr 2>/dev/null)"
        \rm -f $tmp; fail=$((fail+1)); continue
    fi

    g_gtf=$(awk -F'\t' '$3=="gene"' $wd/braker.gtf | wc -l)
    [[ "$g_gtf" == "0" ]] && g_gtf=$(awk -F'\t' '$3=="transcript"' $wd/braker.gtf | wc -l)
    g_g3=$(awk -F'\t' '$3=="gene"' $tmp | wc -l)
    m_g3=$(awk -F'\t' '$3=="mRNA"' $tmp | wc -l)
    n_prot=$(grep -c '^>' $wd/braker.aa)

    if [[ "$m_g3" != "$n_prot" ]]; then
        log_g "  [$sp] FAIL mRNA=$m_g3 != 蛋白数=$n_prot"
        \rm -f $tmp; fail=$((fail+1)); continue
    fi

    mv $tmp $wd/braker.gff3
    made=$((made+1))
    log_g "  [$made] [$sp] OK  gene=$g_g3  mRNA=$m_g3  行数=$(wc -l < $wd/braker.gff3)"
done < $SPECIES_LIST_FILE

log_g ""
log_g "===== 完成 ====="
log_g "  新生成: $made"
log_g "  已存在跳过: $skip"
log_g "  失败: $fail"
log_g "  现有 gff3 总数: $(ls $BATCH_ROOT/*/braker.gff3 2>/dev/null | wc -l) / 127"
