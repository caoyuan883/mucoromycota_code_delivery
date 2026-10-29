import pandas as pd

df = pd.read_csv('sra_detailed_runinfo.csv')
df['AvgLength'] = pd.to_numeric(df['AvgLength'], errors='coerce')

# 唯一有原文支持的过滤条件：Illumina PE (已在查询时限定 Illumina)
# 严格对标硅藻论文：只过滤 PAIRED，其余不做
filtered = df[df['LibraryLayout'] == 'PAIRED'].copy()

# 按 Entrez 返回的自然顺序（即原表顺序），每物种取前 6 个
# 严格对标原文："the top six entries from the Entrez results"
final_runs = filtered.groupby('Species').head(6)

final_runs.to_csv('sci_grade_download_list_v2.csv', index=False)

print("="*50)
print("Scientific Data 级别筛选报告（严格对标文献）")
print("="*50)
print(f"有 RNA-seq 的物种数（BRAKER3 候选）: {final_runs['Species'].nunique()}")
print(f"最终下载清单 SRR 总数: {len(final_runs)}")
print(f"平均每株菌 SRR 数: {len(final_runs)/final_runs['Species'].nunique():.2f}")
print("="*50)
print("注：共培养数据不排除（对标硅藻论文）")
print("注：比对率 <20% 的 library 将在下载后 HISAT2 比对时丢弃")
print("注：读长无最低阈值要求（两篇文献均未设定）")
print("请运行: python filter_sci_grade_v2.py")
