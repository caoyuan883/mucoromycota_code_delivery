#!/bin/bash
# ============================================================================
# BRAKER3 批处理 v3 — Scientific Data 正式版
# ============================================================================
# vs v2 关键差异：
#   ✅ 新目录 11_braker3_batch_v3_scidata/（和 v1/v2 完全隔离）
#   ✅ 物种名后缀 _braker3_v3（避免 augustus_config 冲突）
#   ✅ 新增 --fungus（Mucoromycota 短内含子模型）
#   ✅ 新增 --busco_lineage + --COMPLEASM_PATH（BUSCO 增强，消除批次效应）
#   ✅ 使用 _genomic_softmasked.clean.fa（Step 0 清洗后）
#   ✅ 86 株 BRAKER3 物种表（动态从 species_list_braker3.txt 读取）
#   保持 8 线程 / 4 并行，动态调度（跑完一株立即补下一株）
# ============================================================================

# ---- 环境路径 --------------------------------------------------------------
# [release refactor] site-specific roots are no longer hard-coded. Export them, or
# point BRAKER_SITE_ENV to a file such as config/braker_site.env.example (copied and edited).
# historical: ENV_ROOT=/media/desk16/tly8529/glo_braker_test/08_env
# historical: DATA_ROOT=/media/desk16/tly8529/mucoromycota_annotation_braker
# historical: BUSCO_LINEAGE_DIR=/media/desk16/tly8529/mucoromycota_odb12
if [[ -n "${BRAKER_SITE_ENV:-}" ]]; then
  [[ -f "$BRAKER_SITE_ENV" ]] || { echo "[00_config] BRAKER_SITE_ENV not found: $BRAKER_SITE_ENV" >&2; exit 2; }
  source "$BRAKER_SITE_ENV"
fi
export ENV_ROOT=${ENV_ROOT:?set ENV_ROOT (software environment root)}
export DATA_ROOT=${DATA_ROOT:?set DATA_ROOT (annotation data root)}
export BUSCO_LINEAGE_DIR=${BUSCO_LINEAGE_DIR:?set BUSCO_LINEAGE_DIR (local mucoromycota_odb12 lineage folder)}
export GENOME_DIR=$DATA_ROOT/03_softmasked_genomes
export BAM_DIR=$DATA_ROOT/RNA_SEQ/05_merged_bam
export PROT_DB=$DATA_ROOT/Fungi.clean.fa

# ---- v3 输出目录（和 v1/v2 完全隔离）---------------------------------------
export BATCH_ROOT=$DATA_ROOT/11_braker3_batch_v3_scidata
export LOG_DIR=$BATCH_ROOT/_logs
export STATE_DIR=$BATCH_ROOT/_state
export LOCK_DIR=$BATCH_ROOT/_locks

# ---- BRAKER 版本后缀 -------------------------------------------------------
export SPECIES_SUFFIX=braker3_v3

# ---- 并发资源 --------------------------------------------------------------
# [release refactor] adjustable via the environment or BRAKER_SITE_ENV.
# Example setting (author-confirmed usual runs): 8 threads per job, up to 4 jobs in parallel;
# 1 job was used when resources were limited or for small restarts.
# historical value in this file: THREADS_PER_JOB=8, MAX_PARALLEL=4
export THREADS_PER_JOB=${THREADS_PER_JOB:-8}
export MAX_PARALLEL=${MAX_PARALLEL:-4}

# ---- BRAKER 二进制路径 -----------------------------------------------------
export BRAKER_SCRIPT=$ENV_ROOT/braker2_bundle_transfer/src/BRAKER_src/scripts/braker.pl
export GENEMARK_PATH=$ENV_ROOT/third_party/GeneMark-ETP-main/bin
export BAMTOOLS_PATH=$ENV_ROOT/third_party/bamtools_install/bin
export DIAMOND_PATH=$ENV_ROOT/third_party/diamond/diamond_v2
export PROTHINT_PATH=$ENV_ROOT/third_party/GeneMark-ETP-main/bin/gmes/ProtHint/bin
export TSEBRA_PATH=$ENV_ROOT/braker2_bundle_transfer/src/TSEBRA/bin
export AUGUSTUS_CONFIG_PATH=$ENV_ROOT/augustus_config

# ---- v3 新增：compleasm + BUSCO ---------------------------------------------
export COMPLEASM_PATH=$ENV_ROOT/third_party/compleasm
# BUSCO_LINEAGE_DIR: set above from the site environment

# ---- 86 株物种表（动态从文件读取） ----------------------------------------
# species_list_braker3.txt 由 Step 0 清洗脚本自动生成
SPECIES_LIST_FILE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/species_list_braker3.txt
if [[ -f "$SPECIES_LIST_FILE" ]]; then
  mapfile -t SPECIES_BRAKER3 < "$SPECIES_LIST_FILE"
  export SPECIES_BRAKER3
else
  echo "[00_config] ⚠️  $SPECIES_LIST_FILE 不存在，请先跑 prepare_species_list.sh" >&2
  export SPECIES_BRAKER3=()
fi

# ---- v3 补丁：把 gffread 加进 PATH（best_by_compleasm.py 需要）----
export PATH="$ENV_ROOT/third_party/GeneMark-ETP-main/tools:$PATH"
# ---- v3 补丁 2: miniprot（compleasm 需要）----
export PATH="$ENV_ROOT/third_party/miniprot:$PATH"
# ---- v3 补丁 3: compleasm.py 需要在 PATH（best_by_compleasm 需要）----
export PATH="$ENV_ROOT/third_party/compleasm:$PATH"

# ---- 工具函数 --------------------------------------------------------------
ts() { date '+%F %T'; }
log() { echo "[$(ts)] $*"; }
