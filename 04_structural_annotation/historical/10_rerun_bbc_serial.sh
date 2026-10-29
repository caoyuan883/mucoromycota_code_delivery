#!/bin/bash
# 串行补跑 best_by_compleasm.py
# 用法：nohup bash 10_rerun_bbc_serial.sh > 10_rerun_bbc_serial.log 2>&1 &

set -u

BATCH=/media/desk16/tly8529/mucoromycota_annotation_braker/11_braker3_batch_v3_scidata
ENV_ROOT=/media/desk16/tly8529/glo_braker_test/08_env
ODB=/media/desk16/tly8529/mucoromycota_odb12
COMPLEASM=$ENV_ROOT/third_party/compleasm
export PATH=$COMPLEASM:$PATH
BEST_BY=$ENV_ROOT/braker2_bundle_transfer/src/TSEBRA/bin/best_by_compleasm.py
TSEBRA=$ENV_ROOT/braker2_bundle_transfer/src/TSEBRA/bin/tsebra.py
GET_ANNO=$ENV_ROOT/braker2_bundle_transfer/src/AUGUSTUS_src/scripts/getAnnoFastaFromJoingenes.py
LIST=$BATCH/_rerun_bbc_list.txt
STATE_DIR=$BATCH/_state
LOG_DIR=$BATCH/_rerun_logs
mkdir -p $LOG_DIR

echo "===== 串行补跑 bbc 开始 $(date '+%F %T') ====="
echo "总数: $(wc -l < $LIST) 株"
echo ""

total=$(wc -l < $LIST)
i=0
success=0
failed=0
skipped=0

while IFS= read -r sp; do
    i=$((i + 1))
    [ -z "$sp" ] && continue
    
    workdir="$BATCH/$sp"
    state_file="$STATE_DIR/$sp.state"
    
    # 检查已 DONE_RERUN → 跳过（断点续跑）
    if [ -f "$state_file" ]; then
        cur_state=$(head -1 "$state_file" | awk '{print $1}')
        if [ "$cur_state" = "DONE_RERUN" ] || [ "$cur_state" = "DONE" ]; then
            # DONE 但没 BUSCO_missing 的也跳过（14 株 DONE 里的 10 株）
            if grep -q "BUSCO_missing" "$state_file"; then
                echo "[$i/$total] $sp: 已完成，跳过"
                skipped=$((skipped + 1))
                continue
            fi
        fi
    fi
    
    echo ""
    echo "===== [$i/$total] $sp 开始 $(date '+%F %T') ====="
    
    # 检查产物
    if [ ! -s "$workdir/braker.aa" ] || [ ! -s "$workdir/genome.fa" ] || [ ! -s "$workdir/GeneMark-ETP/genemark.gtf" ]; then
        echo "  ❌ 产物不齐全，跳过"
        echo "FAIL_BBC $(date '+%F %T') reason=missing_inputs" > "$state_file"
        failed=$((failed + 1))
        continue
    fi
    
    # 清理旧的 bbc 目录（避免 tmp 残留干扰）
    if [ -d "$workdir/bbc" ]; then
        \rm -rf "$workdir/bbc"
    fi
    mkdir -p "$workdir/bbc"
    
    # 跑 best_by_compleasm.py，8 线程，30 分钟超时
    start_ts=$(date +%s)
    log_file="$LOG_DIR/${sp}_bbc.log"
    err_file="$LOG_DIR/${sp}_bbc.err"
    
    timeout 1800 /usr/bin/python3 $BEST_BY \
        -m "$workdir/bbc" \
        -d "$workdir" \
        -g "$workdir/genome.fa" \
        -t 8 \
        -p $ODB \
        -y $TSEBRA \
        -f $GET_ANNO \
        > "$log_file" 2> "$err_file"
    rc=$?
    
    end_ts=$(date +%s)
    elapsed=$((end_ts - start_ts))
    elapsed_min=$((elapsed / 60))
    elapsed_sec=$((elapsed % 60))
    
    # 检查产物
    braker_summary="$workdir/bbc/braker/summary.txt"
    
    if [ $rc -eq 0 ] && [ -s "$braker_summary" ]; then
        # 提取 BUSCO_missing (M) 从 summary.txt
        # summary 格式: "S:92.69%, 5487" 等
        busco_missing=$(grep -oE "M:[0-9]+\.[0-9]+%" "$braker_summary" 2>/dev/null | head -1 | sed 's/M://')
        busco_single=$(grep -oE "S:[0-9]+\.[0-9]+%" "$braker_summary" 2>/dev/null | head -1 | sed 's/S://')
        
        # 更新 state
        old_state=$(cat "$state_file" 2>/dev/null | head -1)
        echo "DONE_RERUN $(date '+%F %T') elapsed=${elapsed_min}m${elapsed_sec}s rc=0 BUSCO_missing=$busco_missing BUSCO_single=$busco_single" > "$state_file"
        echo "  ✅ 成功 (${elapsed_min}m${elapsed_sec}s) BUSCO M=$busco_missing S=$busco_single"
        success=$((success + 1))
    else
        # 失败
        old_state=$(cat "$state_file" 2>/dev/null | head -1)
        echo "FAIL_BBC $(date '+%F %T') elapsed=${elapsed_min}m${elapsed_sec}s rc=$rc old_state=$old_state" > "$state_file"
        echo "  ❌ 失败 rc=$rc (${elapsed_min}m${elapsed_sec}s)"
        echo "  err tail:"
        tail -10 "$err_file" 2>/dev/null | sed 's/^/    /'
        failed=$((failed + 1))
    fi
    
    # 每 10 株报告进度
    if [ $((i % 10)) -eq 0 ]; then
        echo ""
        echo "===== 进度 [$i/$total] 成功=$success 失败=$failed 跳过=$skipped ====="
    fi
done < $LIST

echo ""
echo "===== 全部完成 $(date '+%F %T') ====="
echo "总数: $total"
echo "成功: $success"
echo "失败: $failed"
echo "跳过: $skipped"
