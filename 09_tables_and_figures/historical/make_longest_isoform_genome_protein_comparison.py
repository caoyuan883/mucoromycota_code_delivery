#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import re

DESKTOP = Path("/home/py/Desktop")

GENOME_CANDIDATES = [
    DESKTOP / "final_genome_qc_master_table_with_labels.csv",
    (
        DESKTOP
        / "comparison_busco_genome_protein_213"
        / "final_genome_qc_master_table_with_labels.csv"
    ),
]

PROTEIN_CSV = (
    DESKTOP
    / "Mucoromycota_BRAKER_longest_isoform"
    / "longest_isoform_protein_qc_master_table.csv"
)

OUT_DIR = (
    DESKTOP
    / "comparison_busco_genome_longest_isoform_protein_213"
)

WIDE_CSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_comparison_wide.csv"
)

WIDE_TSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_comparison_wide.tsv"
)

LONG_CSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_comparison_long.csv"
)

LONG_TSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_comparison_long.tsv"
)

MU_LONG_CSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_mucoromycota_long.csv"
)

FU_LONG_CSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_fungi_long.csv"
)

DELTA_CSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_delta.csv"
)

REPORT_TXT = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_comparison_report.txt"
)

MATCH_TSV = (
    OUT_DIR
    / "longest_isoform_busco_genome_protein_match_report.tsv"
)

PART_DIR = OUT_DIR / "long_format_4parts"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PART_DIR.mkdir(parents=True, exist_ok=True)

GENOME_CSV = next(
    (
        path
        for path in GENOME_CANDIDATES
        if path.is_file()
    ),
    None,
)

if GENOME_CSV is None:
    raise SystemExit(
        "错误：没有找到final_genome_qc_master_table_with_labels.csv"
    )

if not PROTEIN_CSV.is_file():
    raise SystemExit(
        f"错误：没有找到longest-isoform蛋白QC表：{PROTEIN_CSV}"
    )


def read_csv_file(path):
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames:
            raise ValueError(
                f"{path.name}没有表头"
            )

        headers = [
            field.strip()
            for field in reader.fieldnames
        ]

        rows = []

        for original in reader:
            row = {}

            for old, new in zip(
                reader.fieldnames,
                headers,
            ):
                value = original.get(old, "")
                row[new] = (
                    value.strip()
                    if value is not None
                    else ""
                )

            if any(row.values()):
                rows.append(row)

    return headers, rows


def sample_key(value):
    value = value.strip()

    suffixes = (
        "_genomic_softmasked.fa",
        "_genomic_softmasked.fasta",
        "_longest_isoform.faa",
        "_longest_isoform",
        "_protein.faa",
        ".fasta",
        ".fna",
        ".faa",
        ".fa",
    )

    changed = True

    while changed:
        changed = False

        for suffix in suffixes:
            if value.lower().endswith(
                suffix.lower()
            ):
                value = value[:-len(suffix)]
                changed = True
                break

    return re.sub(
        r"[^a-z0-9]+",
        "",
        value.lower(),
    )


def clean_number(value):
    value = str(value).strip()
    value = value.replace("%", "")
    value = value.replace(",", "")

    if value == "":
        return ""

    try:
        number = float(value)
    except ValueError:
        return ""

    return f"{number:.10g}"


def to_float(value):
    if value in ("", None):
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def find_column(
    headers,
    candidates,
):
    lower = {
        header.lower(): header
        for header in headers
    }

    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]

    return None


def find_busco_column(
    headers,
    metric,
    lineage,
):
    candidates = [
        f"{metric}_{lineage}",
        f"{metric.lower()}_{lineage}",
        f"{lineage}_{metric}",
        f"{lineage}_{metric.lower()}",
    ]

    found = find_column(
        headers,
        candidates,
    )

    if found:
        return found

    normalized_lineage = lineage.lower()
    normalized_metric = metric.lower()

    possible = []

    for header in headers:
        normalized = re.sub(
            r"[^a-z0-9]+",
            "_",
            header.lower(),
        ).strip("_")

        tokens = normalized.split("_")

        if (
            normalized_metric in tokens
            and normalized_lineage in tokens
        ):
            possible.append(header)

    if len(possible) == 1:
        return possible[0]

    return None


genome_headers, genome_rows = read_csv_file(
    GENOME_CSV
)

protein_headers, protein_rows = read_csv_file(
    PROTEIN_CSV
)

genome_sample_col = find_column(
    genome_headers,
    ["sample", "species", "strain", "name"],
)

protein_sample_col = find_column(
    protein_headers,
    ["sample", "species", "strain", "name"],
)

if genome_sample_col is None:
    genome_sample_col = genome_headers[0]

if protein_sample_col is None:
    protein_sample_col = protein_headers[0]

LINEAGES = (
    "mucoromycota",
    "fungi",
)

METRICS = (
    "C",
    "S",
    "D",
    "F",
    "M",
    "n",
)

genome_columns = {}
protein_columns = {}

for lineage in LINEAGES:
    for metric in METRICS:
        genome_columns[
            (lineage, metric)
        ] = find_busco_column(
            genome_headers,
            metric,
            lineage,
        )

        protein_columns[
            (lineage, metric)
        ] = find_busco_column(
            protein_headers,
            metric,
            lineage,
        )

missing_genome_columns = [
    f"{metric}_{lineage}"
    for lineage in LINEAGES
    for metric in METRICS
    if genome_columns[(lineage, metric)] is None
]

missing_protein_columns = [
    f"{metric}_{lineage}"
    for lineage in LINEAGES
    for metric in METRICS
    if protein_columns[(lineage, metric)] is None
]

if missing_genome_columns:
    raise SystemExit(
        "错误：Genome表缺少BUSCO列："
        + ", ".join(missing_genome_columns)
    )

if missing_protein_columns:
    raise SystemExit(
        "错误：Protein表缺少BUSCO列："
        + ", ".join(missing_protein_columns)
    )


def build_index(
    rows,
    sample_column,
    busco_columns,
):
    index = {}
    duplicates = []

    for row in rows:
        sample = row.get(
            sample_column,
            "",
        ).strip()

        if not sample:
            continue

        key_value = sample_key(sample)

        record = {
            "sample": sample,
        }

        for lineage in LINEAGES:
            for metric in METRICS:
                column = busco_columns[
                    (lineage, metric)
                ]

                record[
                    f"{metric}_{lineage}"
                ] = clean_number(
                    row.get(column, "")
                )

        if key_value in index:
            duplicates.append(sample)
            continue

        index[key_value] = record

    return index, duplicates


genome_index, genome_duplicates = build_index(
    genome_rows,
    genome_sample_col,
    genome_columns,
)

protein_index, protein_duplicates = build_index(
    protein_rows,
    protein_sample_col,
    protein_columns,
)

genome_keys = set(genome_index)
protein_keys = set(protein_index)

common_keys = sorted(
    genome_keys & protein_keys,
    key=lambda key_value: (
        genome_index[key_value]["sample"].casefold()
    ),
)

genome_only = sorted(
    genome_keys - protein_keys,
    key=lambda key_value: (
        genome_index[key_value]["sample"].casefold()
    ),
)

protein_only = sorted(
    protein_keys - genome_keys,
    key=lambda key_value: (
        protein_index[key_value]["sample"].casefold()
    ),
)

WIDE_COLUMNS = ["sample"]

for lineage in LINEAGES:
    for dataset in (
        "genome",
        "longest_isoform_protein",
    ):
        for metric in METRICS:
            WIDE_COLUMNS.append(
                f"{metric}_{lineage}_{dataset}"
            )

wide_rows = []
long_rows = []
delta_rows = []
match_rows = []

CATEGORY_NAMES = {
    "S": "Complete (Single-copy)",
    "D": "Complete (Duplicated)",
    "F": "Fragmented",
    "M": "Missing",
}

STACK_ORDER = {
    "S": 1,
    "D": 2,
    "F": 3,
    "M": 4,
}

for sample_order, key_value in enumerate(
    common_keys,
    start=1,
):
    genome = genome_index[key_value]
    protein = protein_index[key_value]

    sample = genome["sample"]

    wide = {
        "sample": sample,
    }

    for lineage in LINEAGES:
        for metric in METRICS:
            wide[
                f"{metric}_{lineage}_genome"
            ] = genome[
                f"{metric}_{lineage}"
            ]

            wide[
                f"{metric}_{lineage}_longest_isoform_protein"
            ] = protein[
                f"{metric}_{lineage}"
            ]

    wide_rows.append(wide)

    for lineage in LINEAGES:
        for dataset, dataset_label, record in (
            (
                "Genome",
                "G",
                genome,
            ),
            (
                "Longest-isoform protein",
                "P",
                protein,
            ),
        ):
            for metric in (
                "S",
                "D",
                "F",
                "M",
            ):
                long_rows.append({
                    "sample": sample,
                    "sample_order": (
                        sample_order
                    ),
                    "lineage": lineage,
                    "dataset": dataset,
                    "dataset_label": (
                        dataset_label
                    ),
                    "metric": metric,
                    "category": (
                        CATEGORY_NAMES[metric]
                    ),
                    "stack_order": (
                        STACK_ORDER[metric]
                    ),
                    "percent": record[
                        f"{metric}_{lineage}"
                    ],
                    "complete_percent": (
                        record[
                            f"C_{lineage}"
                        ]
                    ),
                    "n_buscos": record[
                        f"n_{lineage}"
                    ],
                })

        delta = {
            "sample": sample,
            "lineage": lineage,
        }

        for metric in METRICS:
            genome_value = to_float(
                genome[f"{metric}_{lineage}"]
            )

            protein_value = to_float(
                protein[f"{metric}_{lineage}"]
            )

            delta[
                f"{metric}_genome"
            ] = (
                ""
                if genome_value is None
                else f"{genome_value:.10g}"
            )

            delta[
                f"{metric}_longest_isoform_protein"
            ] = (
                ""
                if protein_value is None
                else f"{protein_value:.10g}"
            )

            if (
                genome_value is not None
                and protein_value is not None
            ):
                delta[
                    f"delta_{metric}_protein_minus_genome"
                ] = f"{protein_value - genome_value:.10g}"
            else:
                delta[
                    f"delta_{metric}_protein_minus_genome"
                ] = ""

        delta_rows.append(delta)

    match_rows.append({
        "sample": sample,
        "genome_sample_original": (
            genome["sample"]
        ),
        "protein_sample_original": (
            protein["sample"]
        ),
        "match_status": "MATCHED",
    })

for key_value in genome_only:
    match_rows.append({
        "sample": (
            genome_index[key_value]["sample"]
        ),
        "genome_sample_original": (
            genome_index[key_value]["sample"]
        ),
        "protein_sample_original": "",
        "match_status": "GENOME_ONLY",
    })

for key_value in protein_only:
    match_rows.append({
        "sample": (
            protein_index[key_value]["sample"]
        ),
        "genome_sample_original": "",
        "protein_sample_original": (
            protein_index[key_value]["sample"]
        ),
        "match_status": "PROTEIN_ONLY",
    })


def write_csv(
    path,
    rows,
    columns,
    delimiter=",",
    bom=True,
):
    encoding = (
        "utf-8-sig"
        if bom
        else "utf-8"
    )

    with path.open(
        "w",
        encoding=encoding,
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=delimiter,
        )
        writer.writeheader()
        writer.writerows(rows)


write_csv(
    WIDE_CSV,
    wide_rows,
    WIDE_COLUMNS,
)

write_csv(
    WIDE_TSV,
    wide_rows,
    WIDE_COLUMNS,
    delimiter="\t",
    bom=False,
)

LONG_COLUMNS = [
    "sample",
    "sample_order",
    "lineage",
    "dataset",
    "dataset_label",
    "metric",
    "category",
    "stack_order",
    "percent",
    "complete_percent",
    "n_buscos",
]

write_csv(
    LONG_CSV,
    long_rows,
    LONG_COLUMNS,
)

write_csv(
    LONG_TSV,
    long_rows,
    LONG_COLUMNS,
    delimiter="\t",
    bom=False,
)

mu_long_rows = [
    row
    for row in long_rows
    if row["lineage"] == "mucoromycota"
]

fu_long_rows = [
    row
    for row in long_rows
    if row["lineage"] == "fungi"
]

write_csv(
    MU_LONG_CSV,
    mu_long_rows,
    LONG_COLUMNS,
)

write_csv(
    FU_LONG_CSV,
    fu_long_rows,
    LONG_COLUMNS,
)

DELTA_COLUMNS = [
    "sample",
    "lineage",
]

for metric in METRICS:
    DELTA_COLUMNS.extend([
        f"{metric}_genome",
        f"{metric}_longest_isoform_protein",
        f"delta_{metric}_protein_minus_genome",
    ])

write_csv(
    DELTA_CSV,
    delta_rows,
    DELTA_COLUMNS,
)

MATCH_COLUMNS = [
    "sample",
    "genome_sample_original",
    "protein_sample_original",
    "match_status",
]

write_csv(
    MATCH_TSV,
    match_rows,
    MATCH_COLUMNS,
    delimiter="\t",
    bom=False,
)

sample_names = [
    genome_index[key_value]["sample"]
    for key_value in common_keys
]

base, remainder = divmod(
    len(sample_names),
    4,
)

part_sizes = [
    base + (
        1
        if index < remainder
        else 0
    )
    for index in range(4)
]

start = 0

for part_number, part_size in enumerate(
    part_sizes,
    start=1,
):
    part_samples = sample_names[
        start:start + part_size
    ]

    part_sample_set = set(part_samples)

    part_rows = [
        row.copy()
        for row in long_rows
        if row["sample"] in part_sample_set
    ]

    local_order = {
        sample: index
        for index, sample in enumerate(
            part_samples,
            start=1,
        )
    }

    for row in part_rows:
        row["sample_order"] = local_order[
            row["sample"]
        ]

    part_file = (
        PART_DIR
        / (
            "longest_isoform_busco_"
            "genome_protein_comparison_"
            f"long_part{part_number}.csv"
        )
    )

    write_csv(
        part_file,
        part_rows,
        LONG_COLUMNS,
    )

    start += part_size

invalid_totals = []

for row in long_rows:
    value = to_float(row["percent"])

    if (
        value is None
        or value < 0
        or value > 100
    ):
        invalid_totals.append(
            (
                row["sample"],
                row["lineage"],
                row["dataset"],
                row["metric"],
                row["percent"],
            )
        )

stack_checks = []

for key_value in common_keys:
    sample = genome_index[
        key_value
    ]["sample"]

    for lineage in LINEAGES:
        for dataset, record in (
            (
                "Genome",
                genome_index[key_value],
            ),
            (
                "Longest-isoform protein",
                protein_index[key_value],
            ),
        ):
            values = [
                to_float(
                    record[
                        f"{metric}_{lineage}"
                    ]
                )
                for metric in (
                    "S",
                    "D",
                    "F",
                    "M",
                )
            ]

            if all(
                value is not None
                for value in values
            ):
                total = sum(values)

                if not (
                    99.7 <= total <= 100.3
                ):
                    stack_checks.append(
                        (
                            sample,
                            lineage,
                            dataset,
                            total,
                        )
                    )

report = [
    "Longest-isoform Genome–Protein BUSCO比较报告\n",
    "============================================\n\n",
    f"Genome表：{GENOME_CSV}\n",
    f"Protein表：{PROTEIN_CSV}\n\n",
    f"Genome唯一样本数：{len(genome_index)}\n",
    f"Protein唯一样本数：{len(protein_index)}\n",
    f"两表共有样本数：{len(common_keys)}\n",
    f"仅Genome表有：{len(genome_only)}\n",
    f"仅Protein表有：{len(protein_only)}\n",
    f"Genome重复sample：{len(genome_duplicates)}\n",
    f"Protein重复sample：{len(protein_duplicates)}\n",
    f"百分比异常记录：{len(invalid_totals)}\n",
    f"S+D+F+M不接近100的组合：{len(stack_checks)}\n",
    "\n仅Genome表有：\n",
]

for key_value in genome_only:
    report.append(
        f"  {genome_index[key_value]['sample']}\n"
    )

report.append("\n仅Protein表有：\n")

for key_value in protein_only:
    report.append(
        f"  {protein_index[key_value]['sample']}\n"
    )

report.append("\nS+D+F+M不接近100：\n")

for (
    sample,
    lineage,
    dataset,
    total,
) in stack_checks:
    report.append(
        f"  {sample}\t{lineage}\t"
        f"{dataset}\t{total:.4f}\n"
    )

REPORT_TXT.write_text(
    "".join(report),
    encoding="utf-8",
)

print("=" * 78)
print("Longest-isoform Genome–Protein BUSCO比较表生成完成")
print("=" * 78)
print(f"Genome表：{GENOME_CSV}")
print(f"Protein表：{PROTEIN_CSV}")
print(f"Genome样本：{len(genome_index)}")
print(f"Protein样本：{len(protein_index)}")
print(f"共有样本：{len(common_keys)}")
print(f"仅Genome：{len(genome_only)}")
print(f"仅Protein：{len(protein_only)}")
print(f"输出目录：{OUT_DIR}")
print(f"宽格式：{WIDE_CSV}")
print(f"总长格式：{LONG_CSV}")
print(f"Mucoromycota长格式：{MU_LONG_CSV}")
print(f"Fungi长格式：{FU_LONG_CSV}")
print(f"变化量表：{DELTA_CSV}")
print(f"四份绘图表：{PART_DIR}")
print(f"报告：{REPORT_TXT}")
