#!/bin/bash
# 02: 诊断卡住的 InterProScan 株（只查看，不停止任何进程）
W=${W:?set W (work directory holding queue.txt, ips/, reports/; historical: /media/desk16/tly8529/func_run)}
IPS=${IPS:?set IPS (path to interproscan.sh; historical: /media/desk16/tly8529/software/interproscan/interproscan-5.78-109.0/interproscan.sh)}
U=${U:-$(id -un)}  # historical default: tly8529
APPL=${APPL:-Pfam,PANTHER,SUPERFAMILY,CDD,SMART,ProSitePatterns,ProSiteProfiles,AntiFam}
NCORE=${NCORE:-72}
mkdir -p $W/reports
TS=$(date '+%Y%m%d_%H%M')
tree(){ ps -u $U -o pid=,ppid= | awk -v r="$1" 'BEGIN{n=split(r,a," ");for(i=1;i<=n;i++)k[a[i]]=1}{p[$1]=$2}END{do{c=0;for(x in p)if(!(x in k)&&(p[x] in k)){k[x]=1;c=1}}while(c);for(x in k)print x}'; }
cmd(){ cat /proc/$1/cmdline 2>/dev/null | tr '\0' ' '; }
pending(){ while read x; do [ -e $W/ips/$x.done ] || echo $x; done < $W/queue.txt; }

cd $W
s=${1:-$(pending | head -1)}
[ -z "$s" ] && { echo "没有未完成的株"; exit 0; }
REP=$W/reports/02_diagnose_${s}_$TS.txt
exec > >(tee $REP) 2>&1
echo "== 02 诊断  $(date '+%F %T')   株: $s"
echo "未完成株列表: $(pending | tr '\n' ' ')"
echo "负载: $(cut -d' ' -f1-3 /proc/loadavg)   核数: $(nproc)"
echo; echo "== 蛋白统计"
echo "蛋白数: $(grep -c '^>' faa/$s.faa)"
echo "最长 10 条（长度  ID）:"
awk '/^>/{if(id)print L"\t"id; id=substr($1,2); L=0; next}{L+=length($0)}END{print L"\t"id}' faa/$s.faa | sort -k1,1nr | head -10
echo "非标准字符（非20种氨基酸）最多的 5 条:"
awk '/^>/{if(id)print x"\t"id; id=substr($1,2); x=0; next}{t=$0; x+=gsub(/[^ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]/,"",t)}END{print x"\t"id}' faa/$s.faa | sort -k1,1nr | head -5
echo; echo "== 进程树（worker 及其子进程）"
for wp in $(pgrep -u $U -f ips_worker.sh); do
  t=$(tree $wp | tr '\n' ',' | sed 's/,$//')
  cmd $(tree $wp | head -1) >/dev/null
  if tree $wp | while read x; do cmd $x; echo; done | grep -q "faa/$s.faa"; then
    echo "-- worker $wp（正在跑 $s）"
    ps -ww -o pid,ppid,etime,time,pcpu,rss,stat,args -p $t | cut -c1-400
    echo "-- 活动子进程中出现的临时文件路径:"
    for x in $(tree $wp); do cmd $x; echo; done | grep -o "[^ ]*tmp/ips_$s[^ ]*" | sort -u | head -20
  fi
done
echo; echo "== 日志末尾 20 行"; tail -20 logs/ips_$s.log
echo; echo "== 日志中的 error/exception/warn（最后 10 条）"; grep -iE "error|exception|warn" logs/ips_$s.log | tail -10
echo; echo "== 临时目录"
du -sh tmp/ips_$s 2>/dev/null
echo "最近 60 分钟内改动过的文件数: $(find tmp/ips_$s -type f -mmin -60 2>/dev/null | wc -l)"
echo "最近 24 小时内改动过的文件数: $(find tmp/ips_$s -type f -mmin -1440 2>/dev/null | wc -l)"
echo "最新改动的 20 个文件:"
find tmp/ips_$s -type f -printf '%TY-%Tm-%Td %TH:%TM  %s  %p\n' 2>/dev/null | sort -r | head -20
echo; echo "== 各成员库在临时目录中的痕迹（按目录名统计文件数）"
find tmp/ips_$s -type f 2>/dev/null | grep -oiE "pfam|panther|superfamily|cdd|smart|prosite[a-z]*|antifam" | sort | uniq -c
echo; echo "== 完成。报告: $REP"
