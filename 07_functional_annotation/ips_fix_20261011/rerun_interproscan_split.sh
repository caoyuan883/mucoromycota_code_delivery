#!/bin/bash
# 03: 停掉卡住的株，切分后并行重跑；超时片段自动二分，定位到单条蛋白后按成员库逐个补跑
W=${W:?set W (work directory holding queue.txt, ips/, reports/; historical: /media/desk16/tly8529/func_run)}
IPS=${IPS:?set IPS (path to interproscan.sh; historical: /media/desk16/tly8529/software/interproscan/interproscan-5.78-109.0/interproscan.sh)}
U=${U:-$(id -un)}  # historical default: tly8529
APPL=${APPL:-Pfam,PANTHER,SUPERFAMILY,CDD,SMART,ProSitePatterns,ProSiteProfiles,AntiFam}
NCORE=${NCORE:-72}
SD=$(cd "$(dirname "$0")" && pwd)
mkdir -p $W/reports
TS=$(date '+%Y%m%d_%H%M')
tree(){ ps -u $U -o pid=,ppid= | awk -v r="$1" 'BEGIN{n=split(r,a," ");for(i=1;i<=n;i++)k[a[i]]=1}{p[$1]=$2}END{do{c=0;for(x in p)if(!(x in k)&&(p[x] in k)){k[x]=1;c=1}}while(c);for(x in k)print x}'; }
cmd(){ cat /proc/$1/cmdline 2>/dev/null | tr '\0' ' '; }
pending(){ while read x; do [ -e $W/ips/$x.done ] || echo $x; done < $W/queue.txt; }

cd $W
s=${1:-$(pending | head -1)}
[ -z "$s" ] && { echo "没有未完成的株"; exit 0; }
if [ -z "$IPS_RERUN_BG" ]; then
  IPS_RERUN_BG=1 setsid nohup bash "$SD/$(basename "$0")" "$s" > $W/reports/03_rerun_${s}_$TS.out 2>&1 < /dev/null &
  echo "已在后台启动，株: $s"
  echo "进度: cat $W/reports/03_progress_$s.txt"
  echo "完整输出: $W/reports/03_rerun_${s}_$TS.out"
  exit 0
fi
R=$W/rerun_$s; PROG=$W/reports/03_progress_$s.txt
log(){ echo "$(date '+%F %T') $*" | tee -a $PROG; }
[ -e ips/$s.done ] && { log "$s 已有 .done，退出"; exit 0; }
[ -d $R ] && { log "!! $R 已存在（之前运行过）。如需重来先手动改名或删除，退出"; exit 1; }
log "== 03 开始  株 $s"

ls $W/reports/02_diagnose_${s}_* >/dev/null 2>&1 || { log "未找到 02 诊断报告，先运行诊断"; IPS_RERUN_BG= bash $SD/diagnose_interproscan.sh $s > /dev/null 2>&1; }

for wp in $(pgrep -u $U -f ips_worker.sh); do
  if tree $wp | while read x; do cmd $x; echo; done | grep -q "faa/$s.faa"; then
    pids=$(tree $wp | tr '\n' ' ')
    log "停止 worker $wp 及其子进程: $pids"
    kill -TERM $pids 2>/dev/null; sleep 30; kill -KILL $pids 2>/dev/null
  fi
done
for cp in $(pgrep -u $U -f ips_ctl.sh); do log "停止控制脚本 $cp"; kill $cp 2>/dev/null; done
sleep 5
[ -e logs/ips_$s.log ] && mv logs/ips_$s.log logs/ips_$s.stuck_$TS.log
\rm -rf tmp/ips_$s
mkdir -p lock/ips_$s

load=$(cut -d' ' -f2 /proc/loadavg)
CCPU=${CCPU:-4}
P=${P:-$(awk -v l=$load -v c=$NCORE -v k=$CCPU 'BEGIN{p=int((c-l-4)/k); if(p<2)p=2; if(p>8)p=8; print p}')}
NCH=${NCH:-16}
export IPS APPL CCPU TMO_PER TMO_MIN
log "负载 $load，并行 $P 个片段，每片段 $CCPU 核，初始切分 $NCH 份"

mkdir -p $R/r1 $R/salvage
cp $SD/run_chunk.sh $R/run_chunk.sh
NSEQ=$(grep -c '^>' faa/$s.faa)
awk -v N=$NCH -v n=$NSEQ -v D=$R/r1 'BEGIN{per=int((n+N-1)/N)} /^>/{k++; c=int((k-1)/per)} {print > sprintf("%s/c%02d_.faa",D,c)}' faa/$s.faa
log "蛋白 $NSEQ 条，切成 $(ls $R/r1/*.faa | wc -l) 份"

: > $R/ok_tsv.list; : > $R/bad_single.list; : > $R/failed_chunks.log
round=1; todo=$(ls $R/r1/*.faa)
while [ -n "$todo" ]; do
  log "第 $round 轮：$(echo "$todo" | wc -l) 个片段"
  echo "$todo" | xargs -P $P -I{} bash $R/run_chunk.sh {}
  next=""; mkdir -p $R/r$((round+1))
  for f in $todo; do
    b=${f%.faa}; rc=99; el=0; T=0; n=$(grep -c '^>' $f); [ -e $b.rc ] && rc=99; el=0; [ -e $b.rc ] && read rc el T n < $b.rc
    if [ "$rc" = 0 ]; then echo $b.tsv >> $R/ok_tsv.list; log "  完成 $(basename $f) n=$n 用时 ${el}s"
    else
      echo "round$round $(basename $f) n=$n rc=$rc elapsed=${el}s timeout=${T}s" >> $R/failed_chunks.log
      log "  失败 $(basename $f) n=$n rc=$rc（124=超时）用时 ${el}s"
      if [ $n -gt 1 ]; then
        bn=$(basename $f .faa); h=$(( (n+1)/2 ))
        awk -v h=$h -v A=$R/r$((round+1))/${bn}a.faa -v B=$R/r$((round+1))/${bn}b.faa '/^>/{k++} {print > (k<=h?A:B)}' $f
        next="$next $R/r$((round+1))/${bn}a.faa $R/r$((round+1))/${bn}b.faa"
      else echo $f >> $R/bad_single.list; fi
    fi
  done
  todo=$(echo $next | tr ' ' '\n' | sed '/^$/d'); round=$((round+1))
done

nb=$(wc -l < $R/bad_single.list)
log "二分结束：无法完整注释的单条蛋白 $nb 条"
: > $R/salvage.tsv.list; printf "protein_id\tanalysis\tstatus\trc\telapsed_s\n" > $R/salvage_summary.tsv
if [ $nb -gt 0 ]; then
  for f in $(cat $R/bad_single.list); do id=$(head -1 $f | cut -c2- | cut -d' ' -f1)
    for a in $(echo $APPL | tr ',' ' '); do echo "$f $R/salvage/$(basename $f .faa).$a $a"; done; done > $R/salvage_jobs.txt
  log "按成员库逐个补跑：$(wc -l < $R/salvage_jobs.txt) 个任务，每个超时 ${TMO_SALV:-3600}s"
  TMO_FIX=${TMO_SALV:-3600} xargs -P $P -n 3 bash $R/run_chunk.sh < $R/salvage_jobs.txt
  while read f b a; do id=$(head -1 $f | cut -c2- | cut -d' ' -f1); rc=99; el=0; [ -e $b.rc ] && read rc el T n < $b.rc
    if [ "$rc" = 0 ]; then echo $b.tsv >> $R/salvage.tsv.list; st=ok; else st=failed; fi
    printf "%s\t%s\t%s\t%s\t%s\n" $id $a $st $rc $el >> $R/salvage_summary.tsv
  done < $R/salvage_jobs.txt
  log "补跑结果："; cat $R/salvage_summary.tsv | tee -a $PROG
fi

nok=0; for b in $(cat $R/ok_tsv.list); do nok=$((nok+$(grep -c '^>' ${b%.tsv}.faa))); done
log "核对：完整完成的蛋白 $nok + 单独处理的蛋白 $nb = $((nok+nb))，原始 $NSEQ"
if [ $((nok+nb)) -ne $NSEQ ]; then log "!! 数量不一致，未生成最终结果，请检查 $R"; exit 1; fi

cat $(cat $R/ok_tsv.list) $(cat $R/salvage.tsv.list) > ips/$s.interproscan.tsv
gzip -f ips/$s.interproscan.tsv
{ echo "株: $s"; echo "原始运行卡住，于 $TS 停止，切分为 $NCH 份重跑（参数与其他株相同: -appl $APPL -goterms -dp -f TSV）"
  echo "完整注释蛋白: $nok / $NSEQ"; echo "需单独处理的蛋白: $nb"
  [ $nb -gt 0 ] && { echo "逐成员库补跑结果:"; cat $R/salvage_summary.tsv; }
  echo "失败片段记录:"; cat $R/failed_chunks.log; } > ips/$s.rerun_note.txt
touch ips/$s.done
log "== 完成：ips/$s.interproscan.tsv.gz  说明: ips/$s.rerun_note.txt"
log "全部株完成数: $(ls ips/*.done | wc -l)/$(wc -l < queue.txt)"
