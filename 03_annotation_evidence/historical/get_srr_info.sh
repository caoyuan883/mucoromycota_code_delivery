#!/bin/bash
INPUT_FILE="sra_availability.csv"
OUTPUT_FILE="sra_detailed_runinfo.csv"

# 写入表头
echo "Species,Run_Accession,LibraryLayout,AvgLength,Size_MB" > $OUTPUT_FILE

echo "正在获取有数据的菌的具体 SRR 编号..."

# 剔除表头和多余的回车符，存入临时文件
tail -n +2 "$INPUT_FILE" | tr -d '\r' > tmp_detailed.csv

# 使用文件描述符 3 读取，避免管道冲突
while IFS=',' read -r species count <&3; do
    if [[ "$count" -gt 0 ]]; then
        # 正常的 esearch 管道传给 efetch
        esearch -db sra -query "\"$species\"[Organism] AND \"illumina\"[Platform] AND \"rna seq\"[Strategy]" | \
        efetch -format runinfo | \
        awk -F',' -v sp="$species" 'NR>1 && $1!="" {print sp","$1","$16","$7","$8}' >> $OUTPUT_FILE
        
        echo "已获取: $species"
        sleep 1
    fi
done 3< "tmp_detailed.csv"

# 清理临时文件
rm tmp_detailed.csv

echo "获取完成！详细的下载编号已保存在 $OUTPUT_FILE"
