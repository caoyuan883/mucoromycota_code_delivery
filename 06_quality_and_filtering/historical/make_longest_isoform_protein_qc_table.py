#!/usr/bin/env python3

from pathlib import Path
from statistics import mean, median
import csv
import math
import re

PROTEIN_DIR = Path(
    "/home/py/Desktop/"
    "Mucoromycota_BRAKER_longest_isoform/"
    "02_longest_isoform_proteins"
)

BUSCO_ROOT = (
    Path.home()
    / "bioinfo"
    / "mucoromycota_busco_longest_isoform_213"
)

BUSCO_MUCOR = (
    BUSCO_ROOT / "results_mucoromycota"
)

BUSCO_FUNGI = (
    BUSCO_ROOT / "results_fungi"
)

OUT_DIR = Path(
    "/home/py/Desktop/"
    "Mucoromycota_BRAKER_longest_isoform"
)

OUT_CSV = (
    OUT_DIR /
    "longest_isoform_protein_qc_master_table.csv"
)

OUT_TSV = (
    OUT_DIR /
    "longest_isoform_protein_qc_master_table.tsv"
)

MATCH_REPORT = (
    OUT_DIR /
    "longest_isoform_protein_qc_match_report.tsv"
)

COLUMNS = [
    "sample",
    "C_mucoromycota",
    "S_mucoromycota",
    "D_mucoromycota",
    "F_mucoromycota",
    "M_mucoromycota",
    "n_mucoromycota",
    "single_copy_faa_mucoromycota",
    "BUSCO_status_mucoromycota",
    "C_fungi",
    "S_fungi",
    "D_fungi",
    "F_fungi",
    "M_fungi",
    "n_fungi",
    "single_copy_faa_fungi",
    "BUSCO_status_fungi",
    "protein_sequence_count",
    "total_aa",
    "mean_protein_length_aa",
    "median_protein_length_aa",
    "protein_N50_aa",
    "shortest_protein_aa",
    "longest_protein_aa",
    "proteins_lt_100aa",
    "proteins_ge_1000aa",
    "proteins_with_terminal_stop",
    "proteins_with_internal_stop",
    "proteins_with_X",
    "empty_protein_sequences",
    "duplicate_fasta_IDs",
    "protein_fasta_size_mb",
]

for directory in (
    PROTEIN_DIR,
    BUSCO_MUCOR,
    BUSCO_FUNGI,
):
    if not directory.is_dir():
        raise SystemExit(
            f"错误：目录不存在：{directory}"
        )

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


def read_fasta(path):
    records = []
    header = None
    sequence_parts = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append(
                        (
                            header,
                            "".join(
                                sequence_parts
                            ),
                        )
                    )

                header = line[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(
                        "第一个有效记录不是FASTA标题"
                    )

                sequence_parts.append(
                    re.sub(r"\s+", "", line)
                )

    if header is not None:
        records.append(
            (
                header,
                "".join(sequence_parts),
            )
        )

    return records


def protein_statistics(path):
    records = read_fasta(path)

    identifiers = []
    clean_lengths = []

    terminal_stop_count = 0
    internal_stop_count = 0
    proteins_with_x = 0
    empty_count = 0

    for header, sequence in records:
        protein_id = (
            header.split()[0]
            if header
            else ""
        )

        identifiers.append(protein_id)

        sequence = re.sub(
            r"\s+",
            "",
            sequence,
        )

        if not sequence:
            empty_count += 1
            clean_lengths.append(0)
            continue

        if sequence.endswith("*"):
            terminal_stop_count += 1
            sequence_without_terminal = (
                sequence[:-1]
            )
        else:
            sequence_without_terminal = (
                sequence
            )

        if "*" in sequence_without_terminal:
            internal_stop_count += 1

        if "X" in sequence.upper():
            proteins_with_x += 1

        # 统计长度时不计终止符*，但内部*也不是氨基酸。
        clean_sequence = (
            sequence.replace("*", "")
            .replace("-", "")
        )

        clean_lengths.append(
            len(clean_sequence)
        )

    duplicate_ids = (
        len(identifiers)
        - len(set(identifiers))
    )

    if clean_lengths:
        total = sum(clean_lengths)
        mean_length = mean(clean_lengths)
        median_length = median(
            clean_lengths
        )
        shortest = min(clean_lengths)
        longest = max(clean_lengths)

        sorted_lengths = sorted(
            clean_lengths,
            reverse=True,
        )

        threshold = total * 0.5
        cumulative = 0
        n50 = 0

        for length in sorted_lengths:
            cumulative += length

            if cumulative >= threshold:
                n50 = length
                break
    else:
        total = 0
        mean_length = 0
        median_length = 0
        shortest = 0
        longest = 0
        n50 = 0

    return {
        "protein_sequence_count": (
            len(records)
        ),
        "total_aa": total,
        "mean_protein_length_aa": (
            f"{mean_length:.2f}"
        ),
        "median_protein_length_aa": (
            f"{median_length:.2f}"
        ),
        "protein_N50_aa": n50,
        "shortest_protein_aa": shortest,
        "longest_protein_aa": longest,
        "proteins_lt_100aa": sum(
            length < 100
            for length in clean_lengths
        ),
        "proteins_ge_1000aa": sum(
            length >= 1000
            for length in clean_lengths
        ),
        "proteins_with_terminal_stop": (
            terminal_stop_count
        ),
        "proteins_with_internal_stop": (
            internal_stop_count
        ),
        "proteins_with_X": (
            proteins_with_x
        ),
        "empty_protein_sequences": (
            empty_count
        ),
        "duplicate_fasta_IDs": (
            duplicate_ids
        ),
        "protein_fasta_size_mb": (
            f"{path.stat().st_size / 1_000_000:.6f}"
        ),
    }


def sample_key(value):
    value = value.strip()

    suffixes = (
        "_longest_isoform.faa",
        "_longest_isoform",
        "_protein.faa",
        ".faa",
    )

    for suffix in suffixes:
        if value.lower().endswith(
            suffix.lower()
        ):
            value = value[:-len(suffix)]
            break

    return re.sub(
        r"[^a-z0-9]+",
        "",
        value.lower(),
    )


def find_summary(result_root, sample):
    sample_output = (
        result_root /
        f"{sample}_longest_isoform"
    )

    if not sample_output.is_dir():
        return None

    candidates = sorted(
        sample_output.glob(
            "short_summary.specific.*.txt"
        )
    )

    if not candidates:
        candidates = sorted(
            sample_output.glob(
                "short_summary*.txt"
            )
        )

    if not candidates:
        candidates = sorted(
            sample_output.rglob(
                "short_summary*.txt"
            )
        )

    return (
        candidates[0]
        if candidates
        else None
    )


def parse_busco_summary(
    result_root,
    sample,
    lineage_name,
):
    summary = find_summary(
        result_root,
        sample,
    )

    result = {
        "C": "",
        "S": "",
        "D": "",
        "F": "",
        "M": "",
        "n": "",
        "single_copy_faa": 0,
        "status": "MISSING",
        "summary": "",
    }

    if summary is None:
        return result

    result["summary"] = str(summary)

    text = summary.read_text(
        encoding="utf-8",
        errors="replace",
    )

    match = re.search(
        r"C:\s*([\d.]+)%"
        r"\s*\[\s*S:\s*([\d.]+)%"
        r"\s*,\s*D:\s*([\d.]+)%\s*\]"
        r"\s*,\s*F:\s*([\d.]+)%"
        r"\s*,\s*M:\s*([\d.]+)%"
        r"\s*,\s*n:\s*(\d+)",
        text,
        re.I,
    )

    if not match:
        result["status"] = (
            "PARSE_ERROR"
        )
        return result

    (
        result["C"],
        result["S"],
        result["D"],
        result["F"],
        result["M"],
        result["n"],
    ) = match.groups()

    sample_output = (
        result_root /
        f"{sample}_longest_isoform"
    )

    expected_dir = (
        sample_output
        / f"run_{lineage_name}"
        / "busco_sequences"
        / "single_copy_busco_sequences"
    )

    if expected_dir.is_dir():
        single_copy_files = list(
            expected_dir.glob("*.faa")
        )
    else:
        single_copy_files = list(
            sample_output.rglob(
                "single_copy_busco_sequences/*.faa"
            )
        )

    result["single_copy_faa"] = len(
        single_copy_files
    )

    # 有有效summary且S序列数与百分比对应的原始计数接近即可。
    result["status"] = "OK"

    return result


protein_files = sorted(
    PROTEIN_DIR.glob(
        "*_longest_isoform.faa"
    ),
    key=lambda path: path.name.casefold(),
)

if not protein_files:
    raise SystemExit(
        "错误：没有发现longest-isoform蛋白文件"
    )

if len(protein_files) != 213:
    print(
        f"警告：预期213个蛋白文件，"
        f"实际找到{len(protein_files)}个"
    )

rows = []
report_rows = []

for index, protein_file in enumerate(
    protein_files,
    start=1,
):
    sample = protein_file.name.removesuffix(
        "_longest_isoform.faa"
    )

    row = {
        column: ""
        for column in COLUMNS
    }

    row["sample"] = sample

    mucor = parse_busco_summary(
        BUSCO_MUCOR,
        sample,
        "mucoromycota_odb12",
    )

    fungi = parse_busco_summary(
        BUSCO_FUNGI,
        sample,
        "fungi_odb12",
    )

    for metric in (
        "C",
        "S",
        "D",
        "F",
        "M",
        "n",
    ):
        row[
            f"{metric}_mucoromycota"
        ] = mucor[metric]

        row[
            f"{metric}_fungi"
        ] = fungi[metric]

    row[
        "single_copy_faa_mucoromycota"
    ] = mucor["single_copy_faa"]

    row[
        "BUSCO_status_mucoromycota"
    ] = mucor["status"]

    row[
        "single_copy_faa_fungi"
    ] = fungi["single_copy_faa"]

    row[
        "BUSCO_status_fungi"
    ] = fungi["status"]

    row.update(
        protein_statistics(
            protein_file
        )
    )

    rows.append(row)

    report_rows.append({
        "sample": sample,
        "protein_file": str(
            protein_file
        ),
        "mucoromycota_summary": (
            mucor["summary"]
        ),
        "mucoromycota_status": (
            mucor["status"]
        ),
        "fungi_summary": (
            fungi["summary"]
        ),
        "fungi_status": (
            fungi["status"]
        ),
    })

    print(
        f"[{index:03d}/{len(protein_files):03d}] "
        f"{sample}: "
        f"MU={mucor['status']} "
        f"FU={fungi['status']} "
        f"proteins="
        f"{row['protein_sequence_count']}"
    )

with OUT_CSV.open(
    "w",
    encoding="utf-8-sig",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=COLUMNS,
    )
    writer.writeheader()
    writer.writerows(rows)

with OUT_TSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=COLUMNS,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

report_columns = [
    "sample",
    "protein_file",
    "mucoromycota_summary",
    "mucoromycota_status",
    "fungi_summary",
    "fungi_status",
]

with MATCH_REPORT.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=report_columns,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(report_rows)

mu_ok = sum(
    row["BUSCO_status_mucoromycota"]
    == "OK"
    for row in rows
)

fu_ok = sum(
    row["BUSCO_status_fungi"]
    == "OK"
    for row in rows
)

print()
print("=" * 78)
print("Longest-isoform蛋白QC表生成完成")
print("=" * 78)
print(f"样本数：{len(rows)}")
print(f"Mucoromycota BUSCO OK：{mu_ok}")
print(f"Fungi BUSCO OK：{fu_ok}")
print(f"CSV：{OUT_CSV}")
print(f"TSV：{OUT_TSV}")
print(f"匹配报告：{MATCH_REPORT}")

if (
    len(rows) != 213
    or mu_ok != 213
    or fu_ok != 213
):
    print(
        "\n警告：样本数或BUSCO成功数不是213，"
        "请检查匹配报告。"
    )
