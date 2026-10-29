#!/bin/bash
# 04: 检查全部 202 株 InterProScan 结果是否完整（只读）
W=${W:-/media/desk16/tly8529/func_run}
U=${U:-tly8529}
cd $W; mkdir -p reports
TS=$(date '+%Y%m%d_%H%M')

CORE="Pfam PANTHER SUPERFAMILY CDD SMART ProSiteProfiles"
OUT=reports/04_ips_check_$TS.tsv; REP=reports/04_ips_check_$TS.txt
printf "strain\tn_proteins\tn_hit\tpct_hit\tn_lines\tanalyses\tmissing_core\tids_not_in_faa\tbad_lines\tgzip_ok\twarn_failed_match\twarn_npe\texit_status\telapsed\tsource\tstatus\n" > $OUT
i=0; N=$(wc -l < queue.txt)
while read s; do i=$((i+1)); printf "\r检查 %d/%d %s          " $i $N $s >&2
  n=$(grep -c '^>' faa/$s.faa); f=ips/$s.interproscan.tsv.gz
  if [ -d rerun_$s ]; then src=rerun; L=$(ls rerun_$s/r*/*.log rerun_$s/salvage/*.log 2>/dev/null)
  else src=normal; L=logs/ips_$s.log; fi
  wf=0; wn=0; for l in $L; do [ -e $l ] || continue; wf=$((wf+$(grep -c "Failed to get match results" $l))); wn=$((wn+$(grep -c "NullPointerException" $l))); done
  if [ $src = normal ]; then ex=$(grep "Exit status:" $L 2>/dev/null | tail -1 | awk '{print $NF}'); el=$(grep "Elapsed (wall clock)" $L 2>/dev/null | tail -1 | awk '{print $NF}')
  else ex=$(cat rerun_$s/r*/*.rc 2>/dev/null | awk '$1!=0{b++}END{print (b?"chunk_fail_"b:"0")}'); el=rerun; fi
  prob=""
  [ -e ips/$s.done ] || prob="$prob;未完成"
  if [ -e $f ] && gzip -t $f 2>/dev/null; then gz=ok
    read nh nl an miss bad <<< $(awk -F'\t' -v core="$CORE" 'NR==FNR{if(/^>/){sub(/^>/,"");split($0,a," ");id[a[1]]=1};next}
      {nl++; if(NF<11)bad++; h[$1]=1; db[$4]=1; if(!($1 in id))ni[$1]=1}
      END{for(k in h)nh++; s="";for(k in db){na++}; m=""; n=split(core,c," "); for(i=1;i<=n;i++) if(!(c[i] in db)) m=m (m?",":"") c[i];
        for(k in ni)nn++; printf "%d %d %d %s %d_%d\n",nh,nl,na,(m?m:"-"),nn+0,bad+0}' faa/$s.faa <(zcat $f))
    nin=${bad%_*}; bl=${bad#*_}
  else gz=bad; nh=0; nl=0; an=0; miss=all; nin=0; bl=0; prob="$prob;gz损坏或缺失"; fi
  [ "$wf" -gt 0 ] && prob="$prob;结果丢失警告$wf"
  [ $src = normal ] && [ "$ex" != 0 ] && prob="$prob;退出码$ex"
  [ $src = rerun ] && [ "$ex" != 0 ] && prob="$prob;重跑片段失败"
  [ "$miss" != "-" ] && prob="$prob;缺少成员库:$miss"
  [ "$nin" -gt 0 ] && prob="$prob;ID不在faa中$nin"
  [ "$bl" -gt 0 ] && prob="$prob;列数异常$bl"
  pct=$(awk -v a=$nh -v b=$n 'BEGIN{printf "%.2f",100*a/b}')
  st=${prob#;}; [ -z "$st" ] && st=OK
  printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" $s $n $nh $pct $nl $an $miss $nin $bl $gz $wf $wn "${ex:-NA}" "${el:-NA}" $src "$st" >> $OUT
done < queue.txt
echo >&2
{ echo "== 04 InterProScan 完整性检查  $(date '+%F %T')"
  echo "株数: $N   .done: $(ls ips/*.done 2>/dev/null | wc -l)"
  awk -F'\t' 'NR>1{n+=$2;h+=$3;if($16=="OK")ok++;else bad++; w+=$11}END{printf "status=OK: %d   有问题: %d\n合计蛋白 %d，有命中 %d（%.2f%%）\n结果丢失警告合计: %d\n",ok,bad,n,h,100*h/n,w}' $OUT
  echo "命中率最低 5 株（仅提示）:"; awk -F'\t' 'NR>1' $OUT | sort -t$'\t' -k4,4n | head -5 | cut -f1-4
  echo "有问题的株:"; awk -F'\t' 'NR>1 && $16!="OK"{print "  "$1"\t"$16}' $OUT
  echo "逐株明细: $W/$OUT"; } | tee $REP
