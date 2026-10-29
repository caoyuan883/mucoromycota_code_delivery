#!/usr/bin/env python3

from pathlib import Path
import csv
import math

ROOT = Path(
    "/home/py/Desktop/comparison_busco_genome_protein_213"
)

INPUT = ROOT / "busco_genome_protein_comparison_long.csv"
OUT_DIR = ROOT / "busco_long_4parts"

if not INPUT.is_file():
    raise SystemExit(f"错误：找不到输入文件：{INPUT}")

OUT_DIR.mkdir(parents=True, exist_ok=True)

with INPUT.open(
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:
    reader = csv.DictReader(handle)

    if not reader.fieldnames:
        raise SystemExit("错误：输入CSV没有表头")

    fieldnames = reader.fieldnames

    if "sample" not in fieldnames:
        raise SystemExit("错误：CSV中没有sample列")

    rows = list(reader)

# 删除sample为空的记录。
rows = [
    row for row in rows
    if row.get("sample", "").strip()
]

# 按sample名称排序，忽略大小写。
rows.sort(
    key=lambda row: (
        row["sample"].strip().casefold(),
        row.get("dataset", "").casefold(),
        int(row.get("stack_order", "0") or 0),
    )
)

# 按排序后的sample获得唯一菌株名单。
samples = []
seen = set()

for row in rows:
    sample = row["sample"].strip()

    if sample not in seen:
        seen.add(sample)
        samples.append(sample)

n_samples = len(samples)

if n_samples == 0:
    raise SystemExit("错误：没有发现有效sample")

# 将样本尽量平均分成4份。
base, remainder = divmod(n_samples, 4)

part_sizes = [
    base + (1 if i < remainder else 0)
    for i in range(4)
]

parts = []
start = 0

for size in part_sizes:
    parts.append(samples[start:start + size])
    start += size

print("=" * 72)
print(f"输入文件：{INPUT}")
print(f"总样本数：{n_samples}")
print(f"总数据行：{len(rows)}")
print(f"每部分样本数：{part_sizes}")
print("=" * 72)

for part_number, part_samples in enumerate(parts, start=1):
    sample_set = set(part_samples)

    part_rows = [
        row for row in rows
        if row["sample"].strip() in sample_set
    ]

    # 重新生成当前部分内部的sample_order。
    order_map = {
        sample: index
        for index, sample in enumerate(part_samples, start=1)
    }

    if "sample_order" in fieldnames:
        for row in part_rows:
            row["sample_order"] = str(
                order_map[row["sample"].strip()]
            )

    output = (
        OUT_DIR /
        f"busco_genome_protein_comparison_long_part{part_number}.csv"
    )

    with output.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(part_rows)

    first_sample = part_samples[0] if part_samples else ""
    last_sample = part_samples[-1] if part_samples else ""

    print(
        f"Part {part_number}: "
        f"{len(part_samples)}个样本，"
        f"{len(part_rows)}行"
    )
    print(f"  范围：{first_sample} → {last_sample}")
    print(f"  文件：{output}")

print("\n完成。每个sample的所有记录均保留在同一部分。")
