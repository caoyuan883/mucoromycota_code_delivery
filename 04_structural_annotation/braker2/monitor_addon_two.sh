#!/bin/bash
# ============================================================================
# 新增 2 株 BRAKER2 精准监控（配合 watch -n 30 使用，只读，绝对安全）
# ============================================================================
# [release refactor] historical: BATCH_ROOT=/media/desk16/tly8529/mucoromycota_annotation_braker/11_braker2_batch_v3_scidata
BATCH_ROOT=${BATCH_ROOT:?set BATCH_ROOT (BRAKER2 batch output folder)}
SPECIES="Cryptocalciella_humilis Rhizopus_sp._TYR1"

echo "==================== $(date '+%F %T') ===================="
echo "磁盘可用: $(df -h "$BATCH_ROOT" | awk 'NR==2{print $4}')   负载:$(uptime | sed 's/.*average://')"
echo ""

for sp in $SPECIES; do
    WD=$BATCH_ROOT/$sp
    ST=$BATCH_ROOT/_state/${sp}.state
    state=$(cat $ST 2>/dev/null || echo "无 state")

    # 已运行时长
    t0=$(echo "$state" | grep -oP '\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}' | head -1)
    if [ -n "$t0" ]; then
        el=$(( $(date +%s) - $(date -d "$t0" +%s) ))
        elapsed=$(printf "%02dh%02dm" $((el/3600)) $(((el%3600)/60)))
    else
        elapsed="--"
    fi

    # 判定当前阶段
    stage="启动中"
    [ -d "$WD/diamond" ]                    && stage="1/6 DIAMOND 蛋白比对"
    [ -d "$WD/Spaln" ]                      && stage="2/6 Spaln/ProtHint 剪接位点"
    [ -s "$WD/prothint.gff" ]               && stage="3/6 ProtHint 完成→GeneMark-EP+"
    [ -s "$WD/GeneMark-EP/genemark.gtf" ]   && stage="4/6 GeneMark 完成→AUGUSTUS 训练"
    [ -s "$WD/augustus.hints.gtf" ]         && stage="5/6 AUGUSTUS 完成→择优/收尾"
    [ -s "$WD/braker.aa" ]                  && stage="6/6 主产物已生成"

    # 进程与 CPU
    pid=$(echo "$state" | grep -oP 'pid=\K\d+')
    if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
        cpu=$(ps -o %cpu= -p $pid 2>/dev/null | tr -d ' ')
        nproc_child=$(pgrep -P $pid 2>/dev/null | wc -l)
        alive="活跃 (主进程 ${cpu}% CPU, 子进程 ${nproc_child})"
    else
        alive="主进程已退出"
    fi
    busy=$(ps -eo pcpu,comm --no-headers | awk '$1>50{print $2}' | sort | uniq -c | sort -rn | head -3 | tr '\n' ' ')

    printf "─── %-26s ───────────────────────\n" "$sp"
    printf "  状态   : %s\n" "$state"
    printf "  已运行 : %s\n" "$elapsed"
    printf "  阶段   : %s\n" "$stage"
    printf "  进程   : %s\n" "$alive"
    printf "  目录   : %s\n" "$(du -sh $WD 2>/dev/null | cut -f1)"

    # 关键产物
    line="  产物   :"
    for f in prothint.gff GeneMark-EP/genemark.gtf augustus.hints.gtf braker.gtf braker.aa; do
        [ -s "$WD/$f" ] && line="$line [$(basename $f)]" 
    done
    [ "$line" = "  产物   :" ] && line="$line （尚无）"
    echo "$line"

    # 最新日志
    last=$(tail -c 4000 "$WD/braker.log" 2>/dev/null | grep -a . | tail -2 | cut -c1-150)
    [ -z "$last" ] && last=$(tail -c 4000 "$BATCH_ROOT/_logs/${sp}.stderr" 2>/dev/null | grep -a . | tail -2 | cut -c1-150)
    echo "  日志   :"
    echo "$last" | sed 's/^/    /'
    echo ""
done

echo "全服 CPU 占用前 3: $(ps -eo pcpu,comm --no-headers | sort -rn | head -3 | awk '{printf "%s(%s%%) ", $2, $1}')"
