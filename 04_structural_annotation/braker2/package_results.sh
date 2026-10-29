#!/bin/bash
# ============================================================================
# Step 07 v4 — 核心产物打包（Scientific Data 投稿版）
# 安全: 只读源目录，只写 _packages_v4，不动旧 _packages
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

MANIFEST=$BATCH_ROOT/_manifest_final.tsv
FULL_ROOT=$BATCH_ROOT/_bbc_full
PKG_DIR=$BATCH_ROOT/_packages_v4
PKG_LOG=$BATCH_ROOT/_package_v4.log
PKG_INDEX=$BATCH_ROOT/_packages_v4_index.tsv

BUSCO_LINEAGE="mucoromycota_odb12"
PIPELINE_VER="BRAKER v2.1.6 (BRAKER2 mode, protein-only evidence)"

mkdir -p $PKG_DIR
log_p() { echo "[$(ts)] [PKG4] $*" | tee -a $PKG_LOG; }

[[ -s "$MANIFEST" ]] || { echo "错误: 找不到 $MANIFEST"; exit 1; }

echo -e "species\tpackage\tsize\tmd5\tn_genes\tn_proteins\tbusco_C\tdecision\tdecision_source" > $PKG_INDEX

log_p "===== 核心产物打包 v4 ====="
log_p "  输出目录: $PKG_DIR"
log_p "  可用空间: $(df -h $BATCH_ROOT | awk 'NR==2{print $4}')"

n=0; ok=0; fail=0
while IFS=$'\t' read -r sp n_gene n_prot C S D F M N bM aM gM dec src nat fset; do
    [[ "$sp" == "species" ]] && continue
    [[ -z "$sp" ]] && continue
    n=$((n+1))

    wd=$BATCH_ROOT/$sp
    vf=$FULL_ROOT/$sp
    stage=$PKG_DIR/$sp

    \rm -rf $stage
    mkdir -p $stage/busco

    miss=""
    for f in braker.gtf braker.gff3 braker.aa braker.codingseq; do
        if [[ -s "$wd/$f" ]]; then cp $wd/$f $stage/; else miss="$miss $f"; fi
    done
    [[ -s "$wd/braker.log" ]] && cp $wd/braker.log $stage/
    [[ -s "$wd/what-to-cite.txt" ]] && cp $wd/what-to-cite.txt $stage/
    [[ -s "$wd/hintsfile.gff" ]] && cp $wd/hintsfile.gff $stage/

    if [[ -n "$miss" ]]; then
        log_p "  [$n/127] [$sp] FAIL 缺文件:$miss"
        fail=$((fail+1)); \rm -rf $stage; continue
    fi

    if [[ "$src" == "braker_native" ]]; then
        cp $wd/best_by_compleasm.log $stage/best_by_compleasm.log
        srcdir=$wd/bbc
        srcnote="Executed automatically by BRAKER during the original run."
    else
        cp $vf.log $stage/best_by_compleasm.log
        srcdir=$vf
        srcnote="BRAKER's internal call failed at runtime (compleasm lineage-file download error); the identical step was re-executed post hoc with the same script, parameters and lineage database."
    fi

    for s in braker augustus genemark; do
        if [[ -s "$srcdir/$s/summary.txt" ]]; then
            mkdir -p $stage/busco/$s
            cp $srcdir/$s/summary.txt $stage/busco/$s/
            [[ -s "$srcdir/$s/full_table.tsv" ]] && cp $srcdir/$s/full_table.tsv $stage/busco/$s/
        fi
    done

    cat > $stage/README.txt <<EOF
================================================================================
 Structural genome annotation: ${sp//_/ }
 Part of: Standardized BRAKER annotations for 127 Mucoromycota genomes
================================================================================

SPECIES              ${sp//_/ }
ANNOTATION PIPELINE  $PIPELINE_VER
EVIDENCE             Protein homology only (no RNA-seq for this genome)
GENERATED            $(date '+%Y-%m-%d')

--------------------------------------------------------------------------------
 GENE SET STATISTICS
--------------------------------------------------------------------------------
 Genes / transcripts        $n_gene
 Predicted proteins         $n_prot
 Final gene set             $fset

--------------------------------------------------------------------------------
 COMPLETENESS  (compleasm, lineage $BUSCO_LINEAGE, N = $N)
--------------------------------------------------------------------------------
 Complete            C = ${C}%
   single-copy       S = ${S}%
   duplicated        D = ${D}%
 Fragmented          F = ${F}%
 Missing             M = ${M}%

--------------------------------------------------------------------------------
 GENE SET SELECTION  (best_by_compleasm)
--------------------------------------------------------------------------------
 Missing BUSCOs per candidate gene set:
   braker            ${bM}%
   augustus          ${aM}%
   genemark          ${gM}%

 Decision            $dec
 Provenance          $src

 $srcnote

 Full command trace and decision statements: best_by_compleasm.log
 Per-set BUSCO tables:                       busco/{braker,augustus,genemark}/

--------------------------------------------------------------------------------
 FILE CONTENTS
--------------------------------------------------------------------------------
 braker.gtf              Final gene structures (GTF)
 braker.gff3             Final gene structures (GFF3)
 braker.aa               Predicted protein sequences (FASTA)
 braker.codingseq        Coding sequences (FASTA)
 hintsfile.gff           Protein-derived extrinsic hints used by BRAKER
 braker.log              Complete BRAKER runtime log
 best_by_compleasm.log   Gene-set selection trace
 busco/                  compleasm summaries and full tables per candidate set
 what-to-cite.txt        Citations required by the BRAKER toolchain
 MD5SUMS                 Checksums of all files in this package

--------------------------------------------------------------------------------
 NOTES
--------------------------------------------------------------------------------
 * Genome assemblies were soft-masked with RepeatModeler/RepeatMasker prior to
   annotation. Assemblies themselves are not redistributed here; accessions are
   listed in the top-level metadata table.
 * All 127 genomes were processed with identical software versions, parameters
   and BUSCO lineage database.
================================================================================
EOF

    cat > $stage/metadata.tsv <<EOF
key	value
species	$sp
pipeline	$PIPELINE_VER
evidence	protein_only
busco_lineage	$BUSCO_LINEAGE
busco_n	$N
n_genes	$n_gene
n_proteins	$n_prot
busco_C	$C
busco_S	$S
busco_D	$D
busco_F	$F
busco_M	$M
braker_missing_pct	$bM
augustus_missing_pct	$aM
genemark_missing_pct	$gM
decision	$dec
decision_source	$src
final_gene_set	$fset
packaged	$(date '+%Y-%m-%d')
EOF

    ( cd $stage && find . -type f ! -name MD5SUMS -print0 | sort -z | xargs -0 md5sum > MD5SUMS )

    tar -czf $PKG_DIR/${sp}_braker2_v4.tar.gz -C $PKG_DIR $sp 2>/dev/null
    rc=$?
    if [[ $rc -ne 0 || ! -s "$PKG_DIR/${sp}_braker2_v4.tar.gz" ]]; then
        log_p "  [$n/127] [$sp] FAIL 压缩失败 rc=$rc"
        fail=$((fail+1)); \rm -rf $stage; continue
    fi
    \rm -rf $stage

    sz=$(du -h $PKG_DIR/${sp}_braker2_v4.tar.gz | cut -f1)
    m5=$(md5sum $PKG_DIR/${sp}_braker2_v4.tar.gz | cut -d' ' -f1)
    echo -e "$sp\t${sp}_braker2_v4.tar.gz\t$sz\t$m5\t$n_gene\t$n_prot\t$C\t$dec\t$src" >> $PKG_INDEX

    ok=$((ok+1))
    log_p "  [$n/127] [$sp] OK  $sz  C=${C}%  $dec"
done < $MANIFEST

log_p ""
log_p "===== 打包完成 ====="
log_p "  成功: $ok    失败: $fail"
log_p "  总大小: $(du -sh $PKG_DIR | cut -f1)"
log_p "  索引:   $PKG_INDEX"
log_p "  剩余空间: $(df -h $BATCH_ROOT | awk 'NR==2{print $4}')"
