#!/bin/bash
INPUT_FILE="12JGI_sra_availability.csv"
OUTPUT_FILE="12JGI_sra_detailed_runinfo.csv"

echo "Species,Run_Accession,LibraryLayout,AvgLength,Size_MB" > "$OUTPUT_FILE"

echo "正在获取有数据的菌的具体 SRR 编号..."

tail -n +2 "$INPUT_FILE" | tr -d '\r' > tmp_12JGI_detailed.csv

while IFS=',' read -r species count <&3; do
  if [[ "$count" -gt 0 ]]; then
    esearch -db sra -query "\"$species\"[Organism] AND \"illumina\"[Platform] AND \"rna seq\"[Strategy]" | \
    efetch -format runinfo | \
    awk -F',' -v sp="$species" 'NR>1 && $1!="" {print sp","$1","$16","$7","$8}' >> "$OUTPUT_FILE"

    echo "已获取: $species"
    sleep 1
  fi
done 3< "tmp_12JGI_detailed.csv"

rm -f tmp_12JGI_detailed.csv

echo "获取完成！详细的下载编号已保存在 $OUTPUT_FILE"
