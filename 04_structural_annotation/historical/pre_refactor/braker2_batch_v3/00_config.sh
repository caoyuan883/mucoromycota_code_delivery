#!/bin/bash
# ============================================================================
# BRAKER2 批处理 v3 — Scientific Data 正式版
# 与 BRAKER3 v3 参数完全对齐：--fungus + --busco_lineage + --COMPLEASM_PATH
# ============================================================================

# ---- 环境路径 --------------------------------------------------------------
export ENV_ROOT=/media/desk16/tly8529/glo_braker_test/08_env
export DATA_ROOT=/media/desk16/tly8529/mucoromycota_annotation_braker

# ---- 输入数据（v3 独立目录，与旧 89 clean 隔离）---------------------------
export GENOME_DIR=$DATA_ROOT/03_softmasked_genomes_v3
export PROT_DB=$DATA_ROOT/Fungi.clean.fa

# ---- v3 输出目录 -----------------------------------------------------------
export BATCH_ROOT=$DATA_ROOT/11_braker2_batch_v3_scidata
export LOG_DIR=$BATCH_ROOT/_logs
export STATE_DIR=$BATCH_ROOT/_state
export LOCK_DIR=$BATCH_ROOT/_locks

# ---- BRAKER 版本后缀（Augustus config 隔离）-------------------------------
export SPECIES_SUFFIX=braker2_v3

# ---- 并发资源（保持 8 线程 × 4 并行）---------------------------------------
export THREADS_PER_JOB=8
export MAX_PARALLEL=1

# ---- BRAKER 及依赖路径 -----------------------------------------------------
export BRAKER_SCRIPT=$ENV_ROOT/braker2_bundle_transfer/src/BRAKER_src/scripts/braker.pl

# GeneMark：BRAKER2 用 gmes_linux_64_4 (EP+ 版 with ProtHint)
export GENEMARK_PATH=$ENV_ROOT/braker2_bundle_transfer/aux/GeneMark/gmes_linux_64_4
export PROTHINT_PATH=$GENEMARK_PATH/ProtHint/bin

# AUGUSTUS
export AUGUSTUS_CONFIG_PATH=$ENV_ROOT/augustus_config
export AUGUSTUS_BIN_PATH=$ENV_ROOT/braker2_bundle_transfer/src/AUGUSTUS_src/bin
export AUGUSTUS_SCRIPTS_PATH=$ENV_ROOT/braker2_bundle_transfer/src/AUGUSTUS_src/scripts

# 第三方工具
export BAMTOOLS_PATH=$ENV_ROOT/third_party/bamtools_install/bin
export DIAMOND_PATH=$ENV_ROOT/third_party/diamond/diamond_v2

# TSEBRA（用了 --busco_lineage 后 BRAKER2 也需要 TSEBRA 做 best_by_compleasm）
# ⚠️ .bashrc 里 TSEBRA_PATH 指向不存在目录，这里显式覆盖为真实路径
export TSEBRA_PATH=$ENV_ROOT/braker2_bundle_transfer/src/TSEBRA/bin

# ---- v3 新增：compleasm + BUSCO（跟 BRAKER3 v3 一致）----------------------
export COMPLEASM_PATH=$ENV_ROOT/third_party/compleasm
export BUSCO_LINEAGE_DIR=/media/desk16/tly8529/mucoromycota_odb12

# PATH：BRAKER 会通过 PATH 找 diamond / miniprot / compleasm
export PATH=$COMPLEASM_PATH:$DIAMOND_PATH:$ENV_ROOT/third_party/miniprot:$PATH

# ---- 125 株 BRAKER2 物种列表 ----------------------------------------------
export SPECIES_LIST_FILE=$BATCH_ROOT/species_list_braker2.txt

# ---- Header 清洗源目录（保留原文件）---------------------------------------
export RAW_GENOME_DIR=$DATA_ROOT/03_softmasked_genomes

# ---- 工具函数 --------------------------------------------------------------
ts() { date '+%F %T'; }
log() { echo "[$(ts)] $*"; }
