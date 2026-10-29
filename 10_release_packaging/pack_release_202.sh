#!/bin/bash
# Pack each strain folder of the release into <Subphylum>/<strain>.tar.gz.
# [release refactor] Paths come from the environment; packing logic unchanged.
#   RELEASE=/path/to/Mucoromycota_202_release SUBPHYLUM_MAP=/path/to/00_subphylum_map.tsv bash pack_release_202.sh
#   (historical MAP location: $R/00_subphylum_map.tsv; derive it with
#    10_release_packaging/derive_subphylum_map.py)
#   STRICT=1 ... -> exit 3 unless 202 archives are OK with 123/57/22 per subphylum
set -uo pipefail
# historical: R=/Users/caoyuan/Desktop/Mucoromycota_202_release
R=${RELEASE:?set RELEASE (folder produced by make_release_202.sh)}
MAP=${SUBPHYLUM_MAP:-$R/00_subphylum_map.tsv}
STRICT=${STRICT:-0}
[ -d "$R" ] || { echo "ERROR: release folder not found: $R" >&2; exit 2; }
[ -s "$MAP" ] || { echo "ERROR: subphylum map not found: $MAP" >&2; exit 2; }
R=$(cd "$R" && pwd)
for sd in Mucoromycotina Mortierellomycotina Glomeromycotina; do
  if [ -d "$R/$sd" ] && [ -n "$(ls -A "$R/$sd" 2>/dev/null)" ] && [ "${OVERWRITE:-0}" != 1 ]; then
    echo "ERROR: $R/$sd already contains archives (set OVERWRITE=1 to repack)" >&2; exit 1
  fi
done
md5list(){ if command -v md5sum >/dev/null 2>&1; then md5sum "$@"; else md5 -r "$@"; fi; }
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
cd "$R" && find Mucoromycotina Mortierellomycotina Glomeromycotina -name "*.tar.gz" | sort | while IFS= read -r f; do md5list "$f"; done > MD5SUMS.txt
echo "压缩包 MD5 已写入 $R/MD5SUMS.txt"
if [ "$STRICT" = 1 ]; then
  awk -F'\t' 'NR>1{c[$2]++; if($5=="OK")ok++; n++}
    END{ if (n!=202 || ok!=202 || c["Mucoromycotina"]!=123 || c["Mortierellomycotina"]!=57 || c["Glomeromycotina"]!=22) exit 3 }' "$LOG" \
    || { echo "STRICT: archive counts/status do not match 202 = 123/57/22 all OK" >&2; exit 3; }
fi
