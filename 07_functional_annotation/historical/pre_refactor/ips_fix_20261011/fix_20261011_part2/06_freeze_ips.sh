#!/bin/bash
# 06: 202 株全部完成且检查通过后，把 InterProScan 结果只读定版到 func_freeze（可加 FORCE=1 跳过警告检查）
SD=$(cd "$(dirname "$0")" && pwd)
W=${W:-/media/desk16/tly8529/func_run}
U=${U:-tly8529}
cd $W; mkdir -p reports
TS=$(date '+%Y%m%d_%H%M')

F=${F:-/media/desk16/tly8529/func_freeze/interproscan_v1_$(date +%Y%m%d)}
N=$(wc -l < queue.txt); nd=$(ls ips/*.done 2>/dev/null | wc -l)
[ $nd -ne $N ] && { echo "!! 只完成 $nd/$N 株，未定版"; exit 1; }
[ -e $F ] && { echo "!! 目标已存在: $F ，未覆盖"; exit 1; }
echo "先运行 04 检查..."; bash $SD/04_check_ips.sh > /dev/null
CK=$(ls -t reports/04_ips_check_*.tsv | head -1); REPK=${CK%.tsv}.txt
nbad=$(awk -F'\t' 'NR>1 && $16!="OK"' $CK | wc -l)
if [ $nbad -gt 0 ] && [ -z "$FORCE" ]; then echo "!! 04 检查有 $nbad 株存在问题，未定版。详见 $REPK"; cat $REPK; exit 1; fi
mkdir -p $F/logs $F/rerun
cp -p ips/*.interproscan.tsv.gz $F/
cp -p ips/*.rerun_note.txt $F/ 2>/dev/null
for s in $(cat queue.txt); do [ -e logs/ips_$s.log ] && cp -p logs/ips_$s.log $F/logs/; done
for d in rerun_*; do [ -d $d ] || continue; mkdir -p $F/rerun/$d
  cp -p $d/*.list $d/*.log $d/*.tsv $F/rerun/$d/ 2>/dev/null
  tar -czf $F/rerun/$d/chunk_logs.tar.gz $d/r*/*.log $d/r*/*.rc $d/salvage/*.log $d/salvage/*.rc 2>/dev/null; done
for l in logs/ips_*.stuck_*.log; do [ -e $l ] && gzip -c $l > $F/logs/$(basename $l).gz; done
cp -p reports/02_diagnose_*.txt reports/03_progress_*.txt $F/rerun/ 2>/dev/null
cp -p $CK $F/ips_check.tsv; cp -p $REPK $F/ips_check_summary.txt
cut -f1-4 $CK > $F/ips_coverage.tsv
cp -p faa.md5 queue.txt $F/
cat > $F/VERSIONS.txt <<'V'
InterProScan 5.78-109.0（本地运行，-dp 关闭预计算查找）
参数: -appl Pfam,PANTHER,SUPERFAMILY,CDD,SMART,ProSitePatterns,ProSiteProfiles,AntiFam -goterms -dp -f TSV
成员库: Pfam 38.2; PANTHER 19.0; SUPERFAMILY 1.75; CDD 3.21; SMART 9.0; PROSITE Patterns/Profiles 2026_01; AntiFam 8.0
Java: OpenJDK 11.0.32.1
输入: 202 株代表性蛋白集（删除全部终止符 *，ID 不变；md5 见 faa.md5）
V
(cd $F && md5sum *.interproscan.tsv.gz ips_check.tsv ips_coverage.tsv faa.md5 queue.txt VERSIONS.txt > MD5SUMS)
chmod -R a-w $F
echo "定版完成: $F"; echo "文件数: $(ls $F | wc -l)   大小: $(du -sh $F | cut -f1)"
cat $F/ips_check_summary.txt
