#!/bin/bash
INPUT_FILE="12JGI_species_list.txt"
OUTPUT_FILE="12JGI_sra_availability.csv"

echo "Species,SRA_Count" > $OUTPUT_FILE
echo "正在自动修复名单换行符并开始查询..."

# 自动处理从 Excel 或 Mac 文本编辑器带来的特殊换行符问题
tr '\r' '\n' < "$INPUT_FILE" | sed '/^$/d' > tmp_list.txt

# 使用文件描述符 3 (<&3) 来读取文件，防止被 NCBI 命令行工具“吞噬”
while IFS= read -r species <&3; do
    # 清除物种名前后多余的空格
    species=$(echo "$species" | xargs)
    if [[ -z "$species" ]]; then continue; fi
    
    # 强制 esearch 不读取外部输入 (< /dev/null)
    count=$(esearch -db sra -query "\"$species\"[Organism] AND \"illumina\"[Platform] AND \"rna seq\"[Strategy]" < /dev/null | xtract -pattern ENTREZ_DIRECT -element Count)
    
    if [[ -z "$count" ]]; then count=0; fi
    
    echo "$species,$count" >> $OUTPUT_FILE
    echo "查询完毕: $species -> $count 个 RNA-seq"
    
    sleep 1
done 3< "tmp_list.txt"

# 删除临时文件
rm tmp_list.txt

echo "完成！全部 211 株查询完毕，结果保存在 $OUTPUT_FILE"
