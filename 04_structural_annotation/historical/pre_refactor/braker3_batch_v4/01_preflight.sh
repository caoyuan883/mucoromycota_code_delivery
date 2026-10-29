#!/bin/bash
# ============================================================================
# Step 1: 起飞前完整检查
# 检查：
#   ✅ 所有 86 株的 clean.fa 存在
#   ✅ 所有 86 株的 merged.bam 存在
#   ✅ 蛋白库 Fungi.clean.fa 可读
#   ✅ 所有 BRAKER 二进制路径可执行
#   ✅ compleasm.py 可运行
#   ✅ BUSCO 库 hmms 目录存在
#   ✅ augustus_config 可写
# ============================================================================
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

LIST=$here/species_list_braker3.txt

echo "════════════════════════════════════════════════"
echo "BRAKER3 v3 preflight"
echo "════════════════════════════════════════════════"

pass=0; fail=0
chk() {
  local msg=$1; local cond=$2
  if eval "$cond" >/dev/null 2>&1; then
    echo "  ✅ $msg"
    pass=$((pass+1))
  else
    echo "  ❌ $msg"
    fail=$((fail+1))
  fi
}

echo ""
echo "-- 输入文件 --"
chk "物种表存在 ($LIST)"                        "[[ -s '$LIST' ]]"
chk "GENOME_DIR 存在 ($GENOME_DIR)"              "[[ -d '$GENOME_DIR' ]]"
chk "BAM_DIR 存在 ($BAM_DIR)"                    "[[ -d '$BAM_DIR' ]]"
chk "PROT_DB 可读 ($PROT_DB)"                    "[[ -s '$PROT_DB' ]]"

echo ""
echo "-- BRAKER 二进制 --"
chk "braker.pl 存在"                             "[[ -x '$BRAKER_SCRIPT' ]]"
chk "GENEMARK gmetp.pl 存在"                     "[[ -f '$GENEMARK_PATH/gmetp.pl' ]]"
chk "PROTHINT prothint.py 存在"                  "[[ -f '$PROTHINT_PATH/prothint.py' ]]"
chk "TSEBRA tsebra.py 存在"                      "[[ -f '$TSEBRA_PATH/tsebra.py' ]]"
chk "DIAMOND 存在"                               "[[ -x '$DIAMOND_PATH/diamond' ]]"
chk "BAMTOOLS 存在"                              "[[ -x '$BAMTOOLS_PATH/bamtools' ]]"
chk "AUGUSTUS_CONFIG_PATH 存在"                  "[[ -d '$AUGUSTUS_CONFIG_PATH' ]]"
chk "AUGUSTUS_CONFIG_PATH 可写"                  "[[ -w '$AUGUSTUS_CONFIG_PATH' ]]"

echo ""
echo "-- v3 新增：compleasm + BUSCO --"
chk "compleasm 目录存在"                          "[[ -d '$COMPLEASM_PATH' ]]"
chk "compleasm.py 存在"                          "[[ -f '$COMPLEASM_PATH/compleasm.py' ]]"
chk "compleasm.py 可运行"                         "python '$COMPLEASM_PATH/compleasm.py' --help >/dev/null 2>&1"
chk "BUSCO 库目录存在"                            "[[ -d '$BUSCO_LINEAGE_DIR' ]]"
chk "BUSCO hmms 目录存在"                         "[[ -d '$BUSCO_LINEAGE_DIR/hmms' ]]"

echo ""
echo "-- 输入文件逐株检查 --"
N_TOTAL=$(wc -l < "$LIST")
miss_fa=0; miss_bam=0
: > /tmp/preflight_missing.txt
while IFS= read -r SP; do
  [[ -z "$SP" ]] && continue
  FA=$GENOME_DIR/${SP}_genomic_softmasked.clean.fa
  BAM=$BAM_DIR/${SP}.merged.bam
  [[ -s "$FA" ]]  || { miss_fa=$((miss_fa+1));  echo "FA_MISS $SP" >> /tmp/preflight_missing.txt; }
  [[ -s "$BAM" ]] || { miss_bam=$((miss_bam+1)); echo "BAM_MISS $SP" >> /tmp/preflight_missing.txt; }
done < "$LIST"
echo "  物种总数: $N_TOTAL"
if [[ $miss_fa -eq 0 && $miss_bam -eq 0 ]]; then
  echo "  ✅ 全部 $N_TOTAL 株的 clean.fa 和 BAM 齐全"
  pass=$((pass+1))
else
  echo "  ❌ 缺失 clean.fa: $miss_fa 株   缺失 BAM: $miss_bam 株"
  echo "     明细：/tmp/preflight_missing.txt"
  cat /tmp/preflight_missing.txt | head -20
  fail=$((fail+1))
fi

echo ""
echo "════════════════════════════════════════════════"
echo "结果: PASS=$pass  FAIL=$fail"
echo "════════════════════════════════════════════════"
if [[ $fail -gt 0 ]]; then
  echo "❌ 请先修复以上 FAIL 项后再启动"
  exit 1
else
  echo "✅ preflight 全部通过，可以启动 supervisor"
  exit 0
fi
