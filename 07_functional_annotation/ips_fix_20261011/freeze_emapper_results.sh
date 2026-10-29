#!/bin/bash
# 01: 检查并定版 eggNOG-mapper 结果（只读复制，不改动原目录）
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
F=${F:-${FREEZE_ROOT:?set FREEZE_ROOT or F (historical FREEZE_ROOT: /media/desk16/tly8529/func_freeze)}/emapper_v1_$(date +%Y%m%d)}
REP=$W/reports/01_freeze_emapper_$TS.txt
exec > >(tee $REP) 2>&1
N=$(wc -l < queue.txt)
echo "== 01 eggNOG 定版  $(date '+%F %T')"
echo "队列株数: $N"
nd=$(ls emapper/*.done 2>/dev/null | wc -l); echo ".done 数: $nd"
nf=$(cat logs/em_failed.txt 2>/dev/null | wc -l); echo "失败记录: $nf"; [ $nf -gt 0 ] && cat logs/em_failed.txt
miss=0; nofoot=0
for s in $(cat queue.txt); do
  [ -s emapper/$s.emapper.annotations ] || { echo "缺少 annotations: $s"; miss=$((miss+1)); continue; }
  [ -e emapper/$s.emapper.seed_orthologs ] || echo "缺少 seed_orthologs: $s"
  tail -5 emapper/$s.emapper.annotations | grep -q "queries scanned" || nofoot=$((nofoot+1))
done
echo "缺少 annotations 的株: $miss"
echo "末尾无 'queries scanned' 统计行的株: $nofoot （若等于 $N，多半只是末尾格式不同，见下方示例）"
echo "-- 示例末尾 3 行:"; tail -3 emapper/$(head -1 queue.txt).emapper.annotations
if [ $nd -ne $N ] || [ $miss -gt 0 ]; then echo "!! 检查未通过，未执行定版"; exit 1; fi
if [ -e $F ]; then echo "!! 目标目录已存在: $F ，未覆盖"; exit 1; fi
mkdir -p $F/logs
cp -p emapper/*.emapper.annotations emapper/*.emapper.seed_orthologs $F/
ls emapper/*.emapper.hits >/dev/null 2>&1 && cp -p emapper/*.emapper.hits $F/
cp -p logs/em_*.log $F/logs/ 2>/dev/null; cp -p logs/emapper_env.txt $F/logs/ 2>/dev/null
cp -p faa.md5 queue.txt $F/
{ printf "strain\tn_proteins\tn_annotated\tpct_annotated\n"
  for s in $(cat queue.txt); do n=$(grep -c '^>' faa/$s.faa); a=$(grep -vc '^#' emapper/$s.emapper.annotations)
    awk -v s=$s -v n=$n -v a=$a 'BEGIN{printf "%s\t%d\t%d\t%.2f\n",s,n,a,100*a/n}'; done; } > $F/emapper_coverage.tsv
awk -F'\t' 'NR>1{n+=$2;a+=$3;p[NR]=$4}END{printf "合计蛋白 %d，有注释 %d（%.2f%%）\n",n,a,100*a/n}' $F/emapper_coverage.tsv
sort -t$'\t' -k4,4n $F/emapper_coverage.tsv | awk 'NR==2{print "最低覆盖: "$0}'
(cd $F && md5sum *.emapper.* emapper_coverage.tsv faa.md5 queue.txt > MD5SUMS)
chmod -R a-w $F
echo "定版目录: $F"; echo "文件数: $(ls $F | wc -l)   大小: $(du -sh $F | cut -f1)"
echo "== 完成。报告: $REP"
