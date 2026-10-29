#!/bin/bash
# ============================================================================
# Step 11 — 顶层 README + 全局 metadata 表（Scientific Data Table 1 素材）
# 只读，只在 _packages_v4/ 下写文档文件
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"

MANIFEST=$BATCH_ROOT/_manifest_final.tsv
PKG_DIR=$BATCH_ROOT/_packages_v4
IDX=$BATCH_ROOT/_packages_v4_index.tsv
T1=$PKG_DIR/Table1_annotation_summary.tsv
T1C=$PKG_DIR/Table1_annotation_summary.csv
RM=$PKG_DIR/README.txt
CKS=$PKG_DIR/MD5SUMS.txt

[[ -s "$MANIFEST" ]] || { echo "缺 $MANIFEST"; exit 1; }
[[ -d "$PKG_DIR" ]]  || { echo "缺 $PKG_DIR"; exit 1; }

echo "[1/4] 生成 Table 1 ..."
echo -e "Species\tAssembly_contig_prefix\tGenes\tProteins\tBUSCO_C_pct\tBUSCO_S_pct\tBUSCO_D_pct\tBUSCO_F_pct\tBUSCO_M_pct\tBUSCO_n\tMissing_braker_pct\tMissing_augustus_pct\tMissing_genemark_pct\tGene_set_selection\tSelection_provenance\tFinal_gene_set\tPackage\tPackage_MD5" > $T1

while IFS=$'\t' read -r sp n_gene n_prot C S D F M N bM aM gM dec src nat fset; do
    [[ "$sp" == "species" ]] && continue
    [[ -z "$sp" ]] && continue
    pfx=$(head -1 $BATCH_ROOT/$sp/braker.gtf 2>/dev/null | cut -f1 | sed 's/[0-9]\{4,\}\.\?[0-9]*$//')
    [[ -z "$pfx" ]] && pfx="NA"
    pkg="${sp}_braker2_v4.tar.gz"
    md5=$(awk -F'\t' -v p="$pkg" '$2==p{print $4}' $IDX 2>/dev/null)
    [[ -z "$md5" ]] && md5="NA"
    echo -e "${sp//_/ }\t$pfx\t$n_gene\t$n_prot\t$C\t$S\t$D\t$F\t$M\t$N\t$bM\t$aM\t$gM\t$dec\t$src\t$fset\t$pkg\t$md5" >> $T1
done < $MANIFEST

python3 - "$T1" "$T1C" <<'PY'
import csv, sys
with open(sys.argv[1]) as f, open(sys.argv[2], 'w', newline='') as o:
    w = csv.writer(o)
    for line in f:
        w.writerow(line.rstrip('\n').split('\t'))
PY

echo "[2/4] 统计 ..."
read NSP MEANC MINC MAXC MEDC TOTP MEANP MINP MAXP TOTG <<< $(awk -F'\t' '
NR>1 {
  n++; c[n]=$5; sc+=$5; sp_+=$4; sg+=$3;
  if(minc==""||$5<minc)minc=$5; if($5>maxc)maxc=$5;
  if(minp==""||$4<minp)minp=$4; if($4>maxp)maxp=$4;
}
END {
  for(i=1;i<=n;i++) for(j=i+1;j<=n;j++) if(c[j]<c[i]){t=c[i];c[i]=c[j];c[j]=t}
  med = (n%2) ? c[int(n/2)+1] : (c[n/2]+c[n/2+1])/2
  printf "%d %.2f %.2f %.2f %.2f %d %.0f %d %d %d", n, sc/n, minc, maxc, med, sp_, sp_/n, minp, maxp, sg
}' $T1)

N90=$(awk -F'\t' 'NR>1 && $5>=90' $T1 | wc -l)
N95=$(awk -F'\t' 'NR>1 && $5>=95' $T1 | wc -l)
N80=$(awk -F'\t' 'NR>1 && $5<80' $T1 | wc -l)
DBEST=$(awk -F'\t' 'NR>1 && $14=="keep_braker_best"' $T1 | wc -l)
DNOADD=$(awk -F'\t' 'NR>1 && $14=="keep_braker_no_add"' $T1 | wc -l)
DMERGE=$(awk -F'\t' 'NR>1 && $14=="merged"' $T1 | wc -l)
SNAT=$(awk -F'\t' 'NR>1 && $15=="braker_native"' $T1 | wc -l)
SPOST=$(awk -F'\t' 'NR>1 && $15=="post_hoc_06c"' $T1 | wc -l)
NPKG=$(ls $PKG_DIR/*.tar.gz 2>/dev/null | wc -l)
TOTSZ=$(du -sh $PKG_DIR 2>/dev/null | cut -f1)

echo "[3/4] 生成顶层 README ..."
cat > $RM <<EOF
================================================================================
 Standardized structural annotations for $NSP Mucoromycota genomes
 generated with BRAKER (protein-evidence mode)
================================================================================

 Release date       $(date '+%Y-%m-%d')
 Genomes annotated  $NSP
 Packages           $NPKG (one gzip-compressed tar archive per genome)
 Total size         $TOTSZ

--------------------------------------------------------------------------------
 1. SUMMARY
--------------------------------------------------------------------------------
This dataset provides structural gene annotations for $NSP fungal genomes of the
phylum Mucoromycota, produced with a single, uniformly parameterized BRAKER
pipeline. All genomes were annotated from protein homology evidence alone; no
RNA-seq data were used for the genomes in this release. Identical software
versions, parameters and BUSCO lineage database were applied throughout, so that
gene sets and completeness metrics are directly comparable across all species.

 Total predicted genes      $TOTG
 Total predicted proteins   $TOTP
 Proteins per genome        mean $MEANP, range $MINP-$MAXP

--------------------------------------------------------------------------------
 2. METHODS
--------------------------------------------------------------------------------
 Repeat masking
   Genome assemblies were soft-masked with RepeatModeler and RepeatMasker before
   annotation.

 Gene prediction
   BRAKER v2.1.6 in BRAKER2 mode (protein-only evidence). Protein hints were
   generated with ProtHint from an OrthoDB-derived Mucoromycota protein set.
   GeneMark-EP+ produced the initial gene models; AUGUSTUS was trained on those
   models and re-predicted genes using the extrinsic protein hints; TSEBRA
   combined the evidence into the final BRAKER gene set.
   Resources per genome: 8 threads, up to 4 genomes annotated in parallel.

 Gene set selection (best_by_compleasm)
   BRAKER's own best_by_compleasm.py compares the completeness of the braker,
   augustus and genemark gene sets and, where the BRAKER set is not the most
   complete one, transfers BUSCO genes that are present only in the AUGUSTUS or
   GeneMark sets into a merged final set (braker_final).

   Outcome across all $NSP genomes:
     keep_braker_best     $DBEST   BRAKER set already the most complete
     keep_braker_no_add   $DNOADD   no additional BUSCOs available to transfer
     merged               $DMERGE   merged set produced
   In this release the BRAKER gene set was retained for every genome; no merged
   gene set was required.

 Completeness assessment
   compleasm (miniprot-based BUSCO implementation), lineage mucoromycota_odb12,
   $(awk -F'\t' 'NR==2{print $10}' $T1) marker genes. The same lineage database copy was used for every
   assessment in this release.

--------------------------------------------------------------------------------
 3. PROVENANCE OF THE GENE-SET SELECTION STEP
--------------------------------------------------------------------------------
 braker_native   $SNAT genomes
   best_by_compleasm was executed automatically by BRAKER during the original
   annotation run.

 post_hoc_06c    $SPOST genomes
   BRAKER's internal call to best_by_compleasm terminated with an error because
   compleasm could not retrieve its lineage placement files (network failure at
   run time). The annotation itself was unaffected. For these genomes the
   identical selection step was re-executed afterwards using the same script,
   the same command-line parameters and the same lineage database, with the
   lineage files provided locally. The provenance of every genome is recorded in
   Table1_annotation_summary.tsv (column Selection_provenance) and in each
   package's metadata.tsv and README.txt.

 GFF3 conversion
   For the $SPOST genomes above, the BRAKER run terminated before the GTF-to-GFF3
   conversion step. braker.gff3 was regenerated with the same AUGUSTUS script
   (gtf2gff.pl --gff3 --printExon --printIntron) used by BRAKER. Equivalence was
   confirmed by re-converting the GTF files of genomes that had completed
   natively and comparing the output to the BRAKER-produced GFF3 (identical).

--------------------------------------------------------------------------------
 4. TECHNICAL VALIDATION
--------------------------------------------------------------------------------
 BUSCO completeness (mucoromycota_odb12)
   mean C     ${MEANC}%
   median C   ${MEDC}%
   range      ${MINC}% - ${MAXC}%
   C >= 95%   $N95 / $NSP genomes
   C >= 90%   $N90 / $NSP genomes
   C <  80%   $N80 / $NSP genomes

 Reproducibility of the completeness assessment
   For the genomes whose selection step had completed natively, the assessment
   was independently repeated with the locally provided lineage database. All
   comparable braker and genemark measurements agreed within 0.5 percentage
   points, confirming that the two execution contexts yield equivalent results.

 Internal consistency
   For every genome the number of mRNA features in braker.gff3 equals the number
   of sequences in braker.aa. Per-file MD5 checksums are included in each
   package (MD5SUMS) and for all archives in MD5SUMS.txt.

--------------------------------------------------------------------------------
 5. FILE STRUCTURE
--------------------------------------------------------------------------------
 README.txt                        this file
 Table1_annotation_summary.tsv     per-genome summary (also provided as .csv)
 MD5SUMS.txt                       checksums of all archives
 <Species>_braker2_v4.tar.gz       one archive per genome, containing:

     braker.gtf              final gene structures (GTF)
     braker.gff3             final gene structures (GFF3)
     braker.aa               predicted proteins (FASTA)
     braker.codingseq        coding sequences (FASTA)
     hintsfile.gff           protein-derived extrinsic hints
     braker.log              complete BRAKER runtime log
     best_by_compleasm.log   gene-set selection trace
     busco/                  compleasm summaries and full tables per candidate set
     README.txt              per-genome documentation
     metadata.tsv            machine-readable per-genome metadata
     what-to-cite.txt        citations required by the BRAKER toolchain
     MD5SUMS                 checksums of the files above

--------------------------------------------------------------------------------
 6. USAGE
--------------------------------------------------------------------------------
   md5sum -c MD5SUMS.txt
   tar -xzf <Species>_braker2_v4.tar.gz
   cd <Species> && md5sum -c MD5SUMS

--------------------------------------------------------------------------------
 7. CITATION
--------------------------------------------------------------------------------
 Please cite the accompanying Data Descriptor together with the tools listed in
 what-to-cite.txt inside each archive (BRAKER, AUGUSTUS, GeneMark-EP+, ProtHint,
 TSEBRA, compleasm/miniprot, RepeatModeler, RepeatMasker).

 Genome assemblies are not redistributed here. Assembly identifiers are given in
 Table1_annotation_summary.tsv (Assembly_contig_prefix).
================================================================================
EOF

echo "[4/4] 生成全局校验和 ..."
( cd $PKG_DIR && md5sum *.tar.gz > MD5SUMS.txt )

echo ""
echo "=========================================="
cat $RM | sed -n '1,12p'
echo "..."
echo "=========================================="
echo ""
echo "Table 1 前 3 行:"
head -3 $T1 | column -t -s$'\t' 2>/dev/null | cut -c1-160
echo ""
echo "生成文件:"
ls -lh $RM $T1 $T1C $CKS | awk '{printf "  %-10s %s\n", $5, $9}'
echo ""
echo "Table 1 行数: $(awk 'NR>1' $T1 | wc -l)   校验和条目: $(wc -l < $CKS)"
