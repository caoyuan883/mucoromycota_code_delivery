#!/bin/bash
# ============================================================================
# Step 07 — 打包核心产物 v2（修复 BUSCO 解析 bug）
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

MODE=${1:-all}

PKG_ROOT=$BATCH_ROOT/_packages
PKG_LOG=$BATCH_ROOT/_package.log
MANIFEST=$BATCH_ROOT/_manifest.tsv

mkdir -p $PKG_ROOT
: > $PKG_LOG

log_pkg() { echo "[$(ts)] [PKG] $*" | tee -a $PKG_LOG; }

log_pkg "===== BRAKER2 v3 核心产物打包 v2 ====="
log_pkg "输出: $PKG_ROOT"
log_pkg "模式: $MODE"

# ---- BUSCO summary.txt 解析函数（正确版）---------------------------------
# 输入: "S:63.91%, 997" → 输出百分比 "63.91"
parse_pct() {
    local line=$1
    # 去掉 "X:" 前缀，去掉 % 及之后所有，去掉前导空格
    echo "$line" | sed 's/^[SDFM]://; s/%.*//; s/^ *//'
}

parse_total() {
    local line=$1
    # "N:1560" → "1560"
    echo "$line" | sed 's/^N://; s/^ *//; s/ *$//'
}

package_one() {
    local sp=$1
    local wd=$BATCH_ROOT/$sp
    local bbc_src=$BATCH_ROOT/_bbc_all/$sp
    local staging=$PKG_ROOT/_staging/$sp
    local tarball=$PKG_ROOT/${sp}_braker2_v3.tar.gz
    
    if [[ ! -s "$wd/braker.aa" ]]; then
        log_pkg "  [$sp] ❌ braker.aa 不存在，跳过"
        echo "$sp SKIP_NO_BRAKER_AA" >> $BATCH_ROOT/_package_skipped.txt
        return 1
    fi
    if [[ ! -s "$wd/braker.gtf" ]]; then
        log_pkg "  [$sp] ❌ braker.gtf 不存在，跳过"
        echo "$sp SKIP_NO_BRAKER_GTF" >> $BATCH_ROOT/_package_skipped.txt
        return 1
    fi
    
    \rm -rf $staging
    mkdir -p $staging/GeneMark-EP $staging/bbc
    
    cp $wd/braker.aa $staging/braker_final.aa
    cp $wd/braker.gtf $staging/braker_final.gtf
    [[ -s "$wd/braker.codingseq" ]] && cp $wd/braker.codingseq $staging/braker_final.codingseq
    
    local has_augustus="no"
    if [[ -s "$wd/augustus.hints.aa" ]]; then
        cp $wd/augustus.hints.aa $staging/
        [[ -s "$wd/augustus.hints.gtf" ]] && cp $wd/augustus.hints.gtf $staging/
        has_augustus="yes"
    fi
    
    if [[ -s "$wd/GeneMark-EP/genemark.gtf" ]]; then
        cp $wd/GeneMark-EP/genemark.gtf $staging/GeneMark-EP/
    fi
    
    [[ -s "$wd/hintsfile.gff" ]] && cp $wd/hintsfile.gff $staging/
    [[ -s "$wd/braker.log" ]] && cp $wd/braker.log $staging/
    [[ -s "$wd/what-to-cite.txt" ]] && cp $wd/what-to-cite.txt $staging/
    
    for sname in braker augustus genemark; do
        local src=$bbc_src/$sname/summary.txt
        if [[ -s "$src" ]]; then
            mkdir -p $staging/bbc/$sname
            cp $src $staging/bbc/$sname/
            [[ -s "$bbc_src/$sname/full_table.tsv" ]] && cp $bbc_src/$sname/full_table.tsv $staging/bbc/$sname/
        fi
    done
    
    # ---- 解析 BUSCO 百分比（braker 集合，用于 manifest/README/state）--------
    local busco_summary=$staging/bbc/braker/summary.txt
    local c_s="NA" c_d="NA" c_f="NA" c_m="NA" c_n="NA"
    if [[ -s "$busco_summary" ]]; then
        local line_s=$(grep '^S:' $busco_summary | head -1)
        local line_d=$(grep '^D:' $busco_summary | head -1)
        local line_f=$(grep '^F:' $busco_summary | head -1)
        local line_m=$(grep '^M:' $busco_summary | head -1)
        local line_n=$(grep '^N:' $busco_summary | head -1)
        
        [[ -n "$line_s" ]] && c_s=$(parse_pct "$line_s")
        [[ -n "$line_d" ]] && c_d=$(parse_pct "$line_d")
        [[ -n "$line_f" ]] && c_f=$(parse_pct "$line_f")
        [[ -n "$line_m" ]] && c_m=$(parse_pct "$line_m")
        [[ -n "$line_n" ]] && c_n=$(parse_total "$line_n")
    fi
    
    local n_prot=$(grep -c '^>' $staging/braker_final.aa)
    
    # ---- 计算 C = S + D（BUSCO 完整度）--------------------------------
    local c_c="NA"
    if [[ "$c_s" != "NA" ]] && [[ "$c_d" != "NA" ]]; then
        c_c=$(echo "$c_s $c_d" | awk '{printf "%.2f", $1+$2}')
    fi
    
    cat > $staging/README.txt <<EOF
Species: $sp
Annotation: BRAKER2 (v3 scidata batch)
Evidence: protein-only (Fungi.clean.fa, 27M proteins)
Used bbc/better (BUSCO-merged) as final: no

Files:
  braker_final.aa       - Final protein set (use this)
  braker_final.gtf      - Final structure
  braker_final.codingseq - Final CDS (if generated)

  braker.log            - BRAKER methodology log
  best_by_compleasm.log - Best gene set decision log (post-hoc)
  hintsfile.gff         - Protein hints (from Fungi.clean.fa)
  augustus.hints.aa/gtf - AUGUSTUS candidate (if available)
  GeneMark-EP/          - GeneMark-EP candidate
  bbc/braker,augustus,genemark/summary.txt - BUSCO reports (mucoromycota_odb12)

BUSCO summary (braker_final.aa):
  Complete (S+D):  ${c_c}%
  Single-copy:     ${c_s}%
  Duplicated:      ${c_d}%
  Fragmented:      ${c_f}%
  Missing:         ${c_m}%
  Total BUSCOs:    $c_n
  Proteins:        $n_prot
EOF
    
    cat > $staging/state.txt <<EOF
DONE $(ts) proteins=$n_prot rc=0 BUSCO_C=${c_c}% BUSCO_S=${c_s}% BUSCO_D=${c_d}% BUSCO_F=${c_f}% BUSCO_M=${c_m}% BUSCO_N=$c_n source_summary=braker has_augustus=$has_augustus
EOF
    
    cat > $staging/best_by_compleasm.log <<EOF
# Post-hoc BUSCO/compleasm evaluation for BRAKER2 v3
# (Reproduces BRAKER3 internal best_by_compleasm behavior for BRAKER2 pipeline)
# Command template:
#   compleasm.py protein -p <aa> -l /media/desk16/tly8529/mucoromycota_odb12 -t 8 -o <out>

EOF
    for sname in genemark braker augustus; do
        local s=$staging/bbc/$sname/summary.txt
        if [[ -s "$s" ]]; then
            local ls=$(grep '^S:' $s | head -1)
            local ld=$(grep '^D:' $s | head -1)
            local lf=$(grep '^F:' $s | head -1)
            local lm=$(grep '^M:' $s | head -1)
            local ss=$(parse_pct "$ls")
            local sd=$(parse_pct "$ld")
            local sf=$(parse_pct "$lf")
            local sm=$(parse_pct "$lm")
            local cc=$(echo "$ss $sd" | awk '{printf "%.2f", $1+$2}')
            echo "[$sname]  Complete=${cc}%  Single=${ss}%  Dup=${sd}%  Frag=${sf}%  Miss=${sm}%" >> $staging/best_by_compleasm.log
        else
            echo "[$sname]  NOT_AVAILABLE (input aa missing)" >> $staging/best_by_compleasm.log
        fi
    done
    echo "" >> $staging/best_by_compleasm.log
    echo "# Decision: keep braker_final.aa (BRAKER2 primary output)" >> $staging/best_by_compleasm.log
    
    cd $PKG_ROOT/_staging
    tar czf $tarball $sp 2>/dev/null
    local rc=$?
    cd - > /dev/null
    
    if [[ "$rc" -ne 0 ]]; then
        log_pkg "  [$sp] ❌ tar 失败 rc=$rc"
        return 1
    fi
    
    \rm -rf $staging
    
    local sz=$(du -h $tarball | awk '{print $1}')
    log_pkg "  [$sp] ✅ $tarball ($sz) prot=$n_prot C=${c_c}%"
    
    echo -e "$sp\t$tarball\t$sz\t$n_prot\t$c_c\t$c_s\t$c_d\t$c_f\t$c_m\t$c_n\t$has_augustus" >> $MANIFEST
    
    return 0
}

echo -e "species\ttarball\tsize\tproteins\tBUSCO_C\tBUSCO_S\tBUSCO_D\tBUSCO_F\tBUSCO_M\tBUSCO_N\thas_augustus" > $MANIFEST
: > $BATCH_ROOT/_package_skipped.txt

i=0
n_total=$(wc -l < $SPECIES_LIST_FILE)

for sp in $(cat $SPECIES_LIST_FILE); do
    i=$((i + 1))
    
    if [[ "$MODE" == "test1" ]] && [[ "$i" -gt 1 ]]; then
        break
    fi
    if [[ "$MODE" == "test5" ]] && [[ "$i" -gt 5 ]]; then
        break
    fi
    
    log_pkg "[$i/$n_total] 打包 $sp"
    package_one "$sp"
done

log_pkg ""
log_pkg "===== 打包完成 ====="
n_pkg=$(ls $PKG_ROOT/*.tar.gz 2>/dev/null | wc -l)
total_sz=$(du -sh $PKG_ROOT | awk '{print $1}')
log_pkg "  已生成 tarball 数: $n_pkg"
log_pkg "  总大小:            $total_sz"
log_pkg "  Manifest:          $MANIFEST"
log_pkg ""
log_pkg "查看 manifest 前 10 行:"
head -10 $MANIFEST | column -t -s $'\t' | tee -a $PKG_LOG

\rm -rf $PKG_ROOT/_staging

if [[ -s $BATCH_ROOT/_package_skipped.txt ]]; then
    log_pkg ""
    log_pkg "⚠️  跳过的物种:"
    cat $BATCH_ROOT/_package_skipped.txt | tee -a $PKG_LOG
fi
