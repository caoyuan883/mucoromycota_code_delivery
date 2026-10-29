#!/bin/bash
# ============================================================================
# Step 02 — BRAKER2 v3 依赖预检
# ============================================================================
set -u
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

fail=0
warn=0

check_file() {
    local desc="$1"
    local path="$2"
    if [ -s "$path" ]; then
        echo "  ✅ $desc"
        echo "     $path"
    else
        echo "  ❌ $desc"
        echo "     $path (缺失或为空)"
        fail=$((fail + 1))
    fi
}

check_dir() {
    local desc="$1"
    local path="$2"
    if [ -d "$path" ]; then
        echo "  ✅ $desc"
        echo "     $path"
    else
        echo "  ❌ $desc"
        echo "     $path (目录不存在)"
        fail=$((fail + 1))
    fi
}

check_exec() {
    local desc="$1"
    local path="$2"
    if [ -x "$path" ]; then
        echo "  ✅ $desc"
        echo "     $path"
    else
        echo "  ❌ $desc"
        echo "     $path (不可执行)"
        fail=$((fail + 1))
    fi
}

echo "===================================="
echo "  BRAKER2 v3 Preflight Check"
echo "  时间: $(ts)"
echo "===================================="

echo ""
echo "===== 1. BRAKER 脚本 ====="
check_file "BRAKER script" "$BRAKER_SCRIPT"
perl_v=$(perl -v 2>/dev/null | head -2 | tail -1)
echo "  Perl: $perl_v"

echo ""
echo "===== 2. GeneMark-EP+ ====="
check_dir "GENEMARK_PATH" "$GENEMARK_PATH"
check_exec "gmes_petap.pl" "$GENEMARK_PATH/gmes_petap.pl"
if [ -f "$GENEMARK_PATH/gm_key_64" ] || [ -f "$HOME/.gm_key" ]; then
    echo "  ✅ GeneMark key 已就位"
else
    echo "  ⚠️  GeneMark key 未找到（在 \$GENEMARK_PATH/gm_key_64 或 \$HOME/.gm_key）"
    warn=$((warn + 1))
fi

echo ""
echo "===== 3. ProtHint ====="
check_dir "PROTHINT_PATH" "$PROTHINT_PATH"
check_exec "prothint.py" "$PROTHINT_PATH/prothint.py"

echo ""
echo "===== 4. AUGUSTUS ====="
check_dir "AUGUSTUS_CONFIG_PATH" "$AUGUSTUS_CONFIG_PATH"
check_dir "AUGUSTUS_BIN_PATH" "$AUGUSTUS_BIN_PATH"
check_dir "AUGUSTUS_SCRIPTS_PATH" "$AUGUSTUS_SCRIPTS_PATH"
check_exec "augustus binary" "$AUGUSTUS_BIN_PATH/augustus"

n_species_config=$(ls $AUGUSTUS_CONFIG_PATH/species/ 2>/dev/null | wc -l)
echo "  当前 augustus_config/species 数量: $n_species_config"

n_v3=$(ls -d $AUGUSTUS_CONFIG_PATH/species/*_braker2_v3 2>/dev/null | wc -l)
if [ "$n_v3" -gt 0 ]; then
    echo "  ⚠️  已存在 $n_v3 个 _braker2_v3 后缀物种（可能是之前测试残留）"
    warn=$((warn + 1))
fi

echo ""
echo "===== 5. compleasm + BUSCO ====="
check_dir "COMPLEASM_PATH" "$COMPLEASM_PATH"
check_exec "compleasm.py" "$COMPLEASM_PATH/compleasm.py"
check_dir "BUSCO_LINEAGE_DIR" "$BUSCO_LINEAGE_DIR"

if [ -d "$BUSCO_LINEAGE_DIR/mucoromycota_odb12" ]; then
    n_hmm=$(ls $BUSCO_LINEAGE_DIR/mucoromycota_odb12/hmms/*.hmm 2>/dev/null | wc -l)
    echo "  ✅ mucoromycota_odb12/hmms: $n_hmm 个 HMM"
elif [ -f "$BUSCO_LINEAGE_DIR/lineages.txt" ] || [ -d "$BUSCO_LINEAGE_DIR/hmms" ]; then
    n_hmm=$(ls $BUSCO_LINEAGE_DIR/hmms/*.hmm 2>/dev/null | wc -l)
    echo "  ✅ BUSCO 库根目录含 HMMs: $n_hmm 个"
else
    echo "  ⚠️  未在 $BUSCO_LINEAGE_DIR 找到标准 BUSCO 结构（hmms/*.hmm）"
    warn=$((warn + 1))
fi

echo ""
echo "===== 6. 第三方工具 ====="
check_dir "BAMTOOLS_PATH" "$BAMTOOLS_PATH"
check_exec "bamtools" "$BAMTOOLS_PATH/bamtools"
check_dir "DIAMOND_PATH" "$DIAMOND_PATH"
check_exec "diamond" "$DIAMOND_PATH/diamond"

echo ""
echo "===== 7. 蛋白库 ====="
check_file "Fungi.clean.fa (OrthoDB v12.2)" "$PROT_DB"
if [ -s "$PROT_DB" ]; then
    n_prot=$(grep -c '^>' $PROT_DB 2>/dev/null || echo "?")
    sz=$(du -h $PROT_DB | awk '{print $1}')
    echo "  蛋白数: $n_prot | 大小: $sz"
fi

echo ""
echo "===== 8. 125 株 clean.fa ====="
n_clean=$(ls $GENOME_DIR/*_genomic_softmasked.clean.fa 2>/dev/null | wc -l)
if [ "$n_clean" -eq 125 ]; then
    echo "  ✅ 125 株 clean.fa 全部就位"
else
    echo "  ❌ clean.fa 数量: $n_clean（期望 125）"
    fail=$((fail + 1))
fi

echo "  抽样验证 header（应无空格）:"
i=0
for fa in $GENOME_DIR/*_genomic_softmasked.clean.fa; do
    i=$((i + 1))
    [ $i -gt 3 ] && break
    header=$(head -1 $fa)
    sp=$(basename $fa _genomic_softmasked.clean.fa)
    if echo "$header" | grep -q ' '; then
        echo "    ❌ $sp: $header (含空格)"
        fail=$((fail + 1))
    else
        echo "    ✅ $sp: $header"
    fi
done

echo ""
echo "===== 9. species 列表 ====="
check_file "species_list_braker2.txt" "$SPECIES_LIST_FILE"
if [ -s "$SPECIES_LIST_FILE" ]; then
    n=$(wc -l < $SPECIES_LIST_FILE)
    echo "  物种数: $n"
fi

echo ""
echo "===== 10. 输出目录准备 ====="
mkdir -p "$BATCH_ROOT" "$LOG_DIR" "$STATE_DIR" "$LOCK_DIR"
check_dir "BATCH_ROOT" "$BATCH_ROOT"
check_dir "LOG_DIR" "$LOG_DIR"
check_dir "STATE_DIR" "$STATE_DIR"
check_dir "LOCK_DIR" "$LOCK_DIR"

echo ""
echo "===== 11. 磁盘空间检查 ====="
df -h $DATA_ROOT | tail -2

echo ""
echo "===================================="
echo "  Preflight 总结"
echo "===================================="
echo "  ❌ 严重错误: $fail"
echo "  ⚠️  警告:     $warn"
echo ""
if [ "$fail" -eq 0 ]; then
    echo "  ✅ 所有依赖就位，可以启动 Step 03 (03_run_braker2_one.sh)"
    if [ "$warn" -gt 0 ]; then
        echo "  ⚠️  但有 $warn 个警告，建议先检查"
    fi
    exit 0
else
    echo "  ❌ 有 $fail 个严重错误，修复后重新运行本预检"
    exit 1
fi
