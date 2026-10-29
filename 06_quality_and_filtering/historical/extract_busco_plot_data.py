#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

root = Path(
    "/home/py/Desktop/comparison_busco_genome_protein_213"
)

genome_file = root / "final_genome_qc_master_table_with_labels.csv"
protein_file = root / "final_protein_qc_master_table.csv"

wide_csv = root / "busco_genome_protein_comparison_wide.csv"
long_csv = root / "busco_genome_protein_comparison_long.csv"
report_file = root / "busco_genome_protein_comparison_report.txt"

for file in (genome_file, protein_file):
    if not file.is_file():
        raise SystemExit(f"文件不存在：{file}")


def clean_sample(value):
    value = str(value).strip()

    suffixes = [
        "_genomic_softmasked.fa",
        "_genomic_softmasked.fasta",
        "_protein.faa",
        ".fasta",
        ".fna",
        ".faa",
        ".fa",
    ]

    for suffix in suffixes:
        if value.lower().endswith(suffix.lower()):
            value = value[:-len(suffix)]
            break

    return value


def make_key(value):
    value = clean_sample(value).lower()
    return re.sub(r"[^a-z0-9]+", "", value)


def clean_percent(series):
    return pd.to_numeric(
        series.astype(str)
        .str.strip()
        .str.replace("%", "", regex=False)
        .str.replace(",", "", regex=False),
        errors="coerce",
    )


genome_raw = pd.read_csv(
    genome_file,
    encoding="utf-8-sig",
)

protein_raw = pd.read_csv(
    protein_file,
    encoding="utf-8-sig",
)

print("Genome原始列：")
for i, column in enumerate(genome_raw.columns, start=1):
    print(f"  {i}: {column}")

print("\nProtein原始列：")
for i, column in enumerate(protein_raw.columns, start=1):
    print(f"  {i}: {column}")

# Genome表：
# 优先按列名提取sample和Mucoromycota C/S/D/F/M。
genome_required = [
    "sample",
    "C_mucoromycota",
    "S_mucoromycota",
    "D_mucoromycota",
    "F_mucoromycota",
    "M_mucoromycota",
]

if all(x in genome_raw.columns for x in genome_required):
    genome = genome_raw[genome_required].copy()
else:
    if genome_raw.shape[1] < 6:
        raise SystemExit(
            "Genome表不足6列，无法提取sample和C/S/D/F/M"
        )

    genome = genome_raw.iloc[:, [0, 1, 2, 3, 4, 5]].copy()
    genome.columns = genome_required

# Protein表：
# 第1列是sample，第4至第8列是C/S/D/F/M。
if protein_raw.shape[1] < 8:
    raise SystemExit(
        "Protein表不足8列，无法提取第1列和第4至第8列"
    )

protein = protein_raw.iloc[:, [0, 3, 4, 5, 6, 7]].copy()

protein.columns = [
    "sample",
    "C_mucoromycota",
    "S_mucoromycota",
    "D_mucoromycota",
    "F_mucoromycota",
    "M_mucoromycota",
]

for table in (genome, protein):
    table["sample"] = table["sample"].map(clean_sample)
    table["key"] = table["sample"].map(make_key)

    for metric in ("C", "S", "D", "F", "M"):
        column = f"{metric}_mucoromycota"
        table[column] = clean_percent(table[column])

genome_duplicates = genome.loc[
    genome["key"].duplicated(keep=False),
    "sample",
].tolist()

protein_duplicates = protein.loc[
    protein["key"].duplicated(keep=False),
    "sample",
].tolist()

genome = genome.drop_duplicates("key", keep="first")
protein = protein.drop_duplicates("key", keep="first")

genome = genome.rename(
    columns={
        "sample": "sample_genome",
        "C_mucoromycota": "C_genome",
        "S_mucoromycota": "S_genome",
        "D_mucoromycota": "D_genome",
        "F_mucoromycota": "F_genome",
        "M_mucoromycota": "M_genome",
    }
)

protein = protein.rename(
    columns={
        "sample": "sample_protein",
        "C_mucoromycota": "C_protein",
        "S_mucoromycota": "S_protein",
        "D_mucoromycota": "D_protein",
        "F_mucoromycota": "F_protein",
        "M_mucoromycota": "M_protein",
    }
)

comparison = genome.merge(
    protein,
    on="key",
    how="outer",
    indicator=True,
)

genome_only = comparison.loc[
    comparison["_merge"] == "left_only",
    "sample_genome",
].dropna().tolist()

protein_only = comparison.loc[
    comparison["_merge"] == "right_only",
    "sample_protein",
].dropna().tolist()

wide = comparison.loc[
    comparison["_merge"] == "both"
].copy()

wide["sample"] = wide["sample_genome"]

wide = wide[
    [
        "sample",
        "C_genome",
        "S_genome",
        "D_genome",
        "F_genome",
        "M_genome",
        "C_protein",
        "S_protein",
        "D_protein",
        "F_protein",
        "M_protein",
    ]
].sort_values(
    "sample",
    key=lambda s: s.str.lower(),
).reset_index(drop=True)

wide.to_csv(
    wide_csv,
    index=False,
    encoding="utf-8-sig",
)

genome_long = wide[
    [
        "sample",
        "C_genome",
        "S_genome",
        "D_genome",
        "F_genome",
        "M_genome",
    ]
].copy()

genome_long.columns = [
    "sample",
    "C",
    "S",
    "D",
    "F",
    "M",
]

genome_long["dataset"] = "Genome"
genome_long["dataset_label"] = "G"

protein_long = wide[
    [
        "sample",
        "C_protein",
        "S_protein",
        "D_protein",
        "F_protein",
        "M_protein",
    ]
].copy()

protein_long.columns = [
    "sample",
    "C",
    "S",
    "D",
    "F",
    "M",
]

protein_long["dataset"] = "Protein"
protein_long["dataset_label"] = "P"

combined = pd.concat(
    [genome_long, protein_long],
    ignore_index=True,
)

long = combined.melt(
    id_vars=[
        "sample",
        "dataset",
        "dataset_label",
        "C",
    ],
    value_vars=["S", "D", "F", "M"],
    var_name="metric",
    value_name="percent",
)

category_names = {
    "S": "Complete (Single-copy)",
    "D": "Complete (Duplicated)",
    "F": "Fragmented",
    "M": "Missing",
}

stack_order = {
    "S": 1,
    "D": 2,
    "F": 3,
    "M": 4,
}

sample_order = {
    sample: i
    for i, sample in enumerate(
        wide["sample"],
        start=1,
    )
}

long["category"] = long["metric"].map(category_names)
long["stack_order"] = long["metric"].map(stack_order)
long["sample_order"] = long["sample"].map(sample_order)

long = long[
    [
        "sample",
        "sample_order",
        "dataset",
        "dataset_label",
        "metric",
        "category",
        "stack_order",
        "percent",
        "C",
    ]
].rename(
    columns={"C": "complete_percent"}
)

long = long.sort_values(
    [
        "sample_order",
        "dataset",
        "stack_order",
    ]
)

long.to_csv(
    long_csv,
    index=False,
    encoding="utf-8-sig",
)

report = [
    "Genome与Protein BUSCO绘图数据提取报告",
    "====================================",
    "",
    f"Genome表样本数：{len(genome)}",
    f"Protein表样本数：{len(protein)}",
    f"共有样本数：{len(wide)}",
    f"仅Genome表有：{len(genome_only)}",
    f"仅Protein表有：{len(protein_only)}",
    "",
    "仅Genome表有：",
]

report.extend(f"  {x}" for x in genome_only)
report.append("")
report.append("仅Protein表有：")
report.extend(f"  {x}" for x in protein_only)

if genome_duplicates:
    report.append("")
    report.append("Genome表重复名称：")
    report.extend(f"  {x}" for x in genome_duplicates)

if protein_duplicates:
    report.append("")
    report.append("Protein表重复名称：")
    report.extend(f"  {x}" for x in protein_duplicates)

report_file.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8",
)

print("\n完成")
print(f"共有样本：{len(wide)}")
print(f"仅Genome：{len(genome_only)}")
print(f"仅Protein：{len(protein_only)}")
print(f"宽格式：{wide_csv}")
print(f"长格式：{long_csv}")
print(f"报告：{report_file}")
