#!/bin/bash
# ============================================================================
# 单株 BRAKER2 精准监控（含 CPU / 内存）
# 用法: watch -n 30 -t 'bash ~/braker2_batch_v3/18_watch_one.sh Umbelopsis_ramanniana'
# 只读，绝对安全
# ============================================================================
BR=/media/desk16/tly8529/mucoromycota_annotation_braker/11_braker2_batch_v3_scidata
SP=${1:-Umbelopsis_ramanniana}
WD=$BR/$SP
ST=$BR/_state/${SP}.state

echo "══════════ $SP ══════════  $(date '+%F %T')"

state=$(cat $ST 2>/dev/null || echo "无 state")
echo "状态: $state"

# 已运行时长
t0=$(echo "$state" | grep -oP '\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}' | head -1)
if [ -n "$t0" ]; then
    el=$(( $(date +%s) - $(date -d "$t0" +%s) ))
    printf "已运行: %02dh%02dm%02ds\n" $((el/3600)) $(((el%3600)/60)) $((el%60))
fi

# 阶段判定
stage="1/6 DIAMOND 蛋白比对"
[ -d "$WD/Spaln" ]                    && stage="2/6 Spaln / ProtHint 剪接位点"
[ -s "$WD/prothint.gff" ]             && stage="3/6 ProtHint 完成 → GeneMark-EP+"
[ -s "$WD/GeneMark-EP/genemark.gtf" ] && stage="4/6 GeneMark 完成 → AUGUSTUS 训练"
[ -s "$WD/augustus.hints.gtf" ]       && stage="5/6 AUGUSTUS 完成 → 择优收尾"
[ -s "$WD/braker.aa" ]                && stage="6/6 主产物已生成"
echo "阶段: $stage"
echo "目录: $(du -sh $WD 2>/dev/null | cut -f1)"

# 产物
line="产物:"
for f in prothint.gff GeneMark-EP/genemark.gtf augustus.hints.gtf braker.gtf braker.aa braker.gff3; do
    [ -s "$WD/$f" ] && line="$line [$(basename $f)]"
done
[ "$line" = "产物:" ] && line="$line （尚无）"
echo "$line"

# ---------------- 本任务的进程树 CPU / 内存 ----------------
pid=$(echo "$state" | grep -oP 'pid=\K\d+')
echo ""
echo "────── 本任务资源占用 ──────"
if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
    tree=$(pstree -p $pid 2>/dev/null | grep -oP '\(\K\d+' | tr '\n' ',' | sed 's/,$//')
    [ -z "$tree" ] && tree=$pid
    read tcpu tmem trss nproc <<<$(ps -o pcpu=,pmem=,rss= -p "$tree" 2>/dev/null | \
        awk '{c+=$1; m+=$2; r+=$3; n++} END{printf "%.1f %.1f %.1f %d", c, m, r/1048576, n}')
    echo "主进程 pid=$pid  进程数=$nproc"
    echo "CPU 合计: ${tcpu:-0}%   内存合计: ${tmem:-0}%  (${trss:-0} GB RSS)"
    echo "占用最高的 5 个子进程:"
    ps -o pcpu=,pmem=,rss=,etime=,comm= -p "$tree" 2>/dev/null | sort -rn | head -5 | \
        awk '{printf "   %6s%%CPU %5s%%MEM %6.2fGB %10s  %s\n", $1, $2, $3/1048576, $4, $5}'
else
    echo "主进程已退出（若状态仍为 RUNNING 则异常）"
    echo "残留相关进程: $(pgrep -fc "$SP" 2>/dev/null || echo 0) 个"
fi

# ---------------- 全机状况 ----------------
echo ""
echo "────── 服务器整体 ──────"
ncpu=$(nproc)
echo "负载: $(cut -d' ' -f1-3 /proc/loadavg)   逻辑核: $ncpu"
free -g | awk 'NR==2{printf "内存: 已用 %d GB / 共 %d GB (%.0f%%)   可用 %d GB\n", $3, $2, 100*$3/$2, $7}'
echo "磁盘: /media/desk16 可用 $(df -h /media/desk16 | awk 'NR==2{print $4}')"
echo "全机 CPU 前 5:"
ps -eo pcpu,pmem,user,comm --no-headers | sort -rn | head -5 | \
    awk '{printf "   %6s%%CPU %5s%%MEM  %-12s %s\n", $1, $2, $3, $4}'

# ---------------- 日志 ----------------
echo ""
echo "────── braker.log 末尾 ──────"
tail -c 8000 "$WD/braker.log" 2>/dev/null | grep -a . | tail -4 | cut -c1-135 | sed 's/^/  /'
err=$(tail -c 4000 "$BR/_logs/${SP}.stderr" 2>/dev/null | grep -aiE 'error|failed|cannot|died' | tail -2 | cut -c1-135)
[ -n "$err" ] && { echo "  ⚠ stderr 异常:"; echo "$err" | sed 's/^/    /'; }
