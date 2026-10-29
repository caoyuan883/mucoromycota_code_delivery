#!/bin/bash
set -uo pipefail
R=/Users/caoyuan/Desktop/Mucoromycota_202_release
MAP=$R/00_subphylum_map.tsv
export COPYFILE_DISABLE=1
mkdir -p "$R/Mucoromycotina" "$R/Mortierellomycotina" "$R/Glomeromycotina"
LOG=$R/00_pack_log.tsv
printf "strain\tsubphylum\tarchive\tn_files\tstatus\n" > "$LOG"

for d in "$R"/*/; do
  sp=$(basename "$d")
  case "$sp" in Mucoromycotina|Mortierellomycotina|Glomeromycotina) continue;; esac
  code=$(awk -F'\t' -v s="$sp" '$1==s{print $2}' "$MAP")
  case "$code" in Muc) sub=Mucoromycotina;; Mor) sub=Mortierellomycotina;; Glo) sub=Glomeromycotina;;
    *) printf "%s\t%s\t-\t-\tNO_SUBPHYLUM\n" "$sp" "$code" >> "$LOG"; continue;; esac
  ok=1
  for f in _longest_isoform.gtf _longest_isoform.gff3 _longest_isoform.faa _interproscan.tsv _emapper.annotations.tsv; do
    [ -s "$d/${sp}$f" ] || ok=0
  done
  [ $ok -eq 1 ] || { printf "%s\t%s\t-\t-\tMISSING_FILE\n" "$sp" "$sub" >> "$LOG"; continue; }
  a=$R/$sub/${sp}.tar.gz
  tar --exclude '.DS_Store' -czf "$a" -C "$R" "$sp"
  n=$(tar -tzf "$a" | grep -v '/$' | wc -l | tr -d ' ')
  if gzip -t "$a" && [ "$n" -eq 5 ]; then st=OK; else st=FAIL; fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$sp" "$sub" "$sub/${sp}.tar.gz" "$n" "$st" >> "$LOG"
  echo "$st  $sub/$sp"
done

echo "== 汇总 =="
awk -F'\t' 'NR>1{c[$2]++; s[$5]++} END{for(k in c)print "  "k": "c[k]; for(k in s)print "  状态 "k": "s[k]}' "$LOG"
echo "期望: Mucoromycotina 123, Mortierellomycotina 57, Glomeromycotina 22；状态全为 OK 202"
cd "$R" && find Mucoromycotina Mortierellomycotina Glomeromycotina -name "*.tar.gz" | sort | xargs md5 -r > MD5SUMS.txt
echo "压缩包 MD5 已写入 $R/MD5SUMS.txt"
