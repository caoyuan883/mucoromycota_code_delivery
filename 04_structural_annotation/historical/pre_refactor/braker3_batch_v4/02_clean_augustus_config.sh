#!/bin/bash
# ============================================================================
# Step 2: 清理 augustus_config/species/ 里的残留物种目录
# 只清理带 _braker3_v3 后缀的目录，绝不动其它任何东西
#
# ⚠️ 关键：
#   - 使用 \rm 或 /bin/rm 绕过 trash-put 别名
#   - 只匹配 *_braker3_v3 后缀，不动其它
#   - 手动确认后才执行
# ============================================================================
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/00_config.sh"

SPECIES_ROOT=$AUGUSTUS_CONFIG_PATH/species

if [[ ! -d "$SPECIES_ROOT" ]]; then
  echo "❌ augustus species 目录不存在: $SPECIES_ROOT"
  exit 1
fi

echo "════════════════════════════════════════════════"
echo "清理 augustus_config/species/*_braker3_v3"
echo "════════════════════════════════════════════════"
echo "目录: $SPECIES_ROOT"
echo ""

# 列出待清理目录
MATCHES=()
for d in "$SPECIES_ROOT"/*_${SPECIES_SUFFIX}; do
  [[ -d "$d" ]] || continue
  MATCHES+=("$d")
done

if [[ ${#MATCHES[@]} -eq 0 ]]; then
  echo "✅ 无需清理（没有 *_${SPECIES_SUFFIX} 残留目录）"
  exit 0
fi

echo "发现 ${#MATCHES[@]} 个残留目录："
for d in "${MATCHES[@]}"; do
  echo "  - $(basename "$d")"
done
echo ""

read -p "确认清理？(输入 yes 继续) " ans
if [[ "$ans" != "yes" ]]; then
  echo "已取消"
  exit 0
fi

# 使用 \rm 绕过 trash-put 别名
for d in "${MATCHES[@]}"; do
  echo "清理: $d"
  \rm -rf "$d"
done

echo "✅ 清理完成"
