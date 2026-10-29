#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv
import re
import shutil
import sys

RESULT_ROOT = Path(
    "/home/py/bioinfo/"
    "mucoromycota_busco_longest_isoform_213/"
    "results_mucoromycota"
)

OUTPUT_ROOT = Path(
    "/home/py/Desktop/"
    "Mucoromycota_213_single_copy_BUSCO"
)

BY_SPECIES = OUTPUT_ROOT / "by_species"
BY_BUSCO = OUTPUT_ROOT / "by_BUSCO"

REPORT = OUTPUT_ROOT / "extraction_report.tsv"
PRESENCE = OUTPUT_ROOT / "BUSCO_presence_counts.tsv"
SAMPLE_LIST = OUTPUT_ROOT / "sample_list.txt"

for directory in (
    BY_SPECIES,
    BY_BUSCO,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

if not RESULT_ROOT.is_dir():
    raise SystemExit(
        f"错误：BUSCO结果目录不存在：{RESULT_ROOT}"
    )


def read_fasta(path):
    records = []
    header = None
    sequence = []

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
                            "".join(sequence),
                        )
                    )

                header = line[1:].strip()
                sequence = []
            else:
                if header is None:
                    raise ValueError(
                        "FASTA首个有效记录不是标题"
                    )

                sequence.append(
                    re.sub(r"\s+", "", line)
                )

    if header is not None:
        records.append(
            (
                header,
                "".join(sequence),
            )
        )

    return records


def wrap(sequence, width=60):
    return "\n".join(
        sequence[index:index + width]
        for index in range(
            0,
            len(sequence),
            width,
        )
    )


def write_species_fasta(
    path,
    sample,
    busco_id,
    records,
):
    """
    每个输出FASTA标题以物种名为主。
    正常单拷贝文件应只有一条记录。
    """
    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for index, (
            original_header,
            sequence,
        ) in enumerate(
            records,
            start=1,
        ):
            if len(records) == 1:
                new_header = (
                    f"{sample} BUSCO={busco_id}"
                )
            else:
                new_header = (
                    f"{sample}__copy{index} "
                    f"BUSCO={busco_id}"
                )

            handle.write(
                f">{new_header}\n"
            )
            handle.write(
                wrap(sequence) + "\n"
            )


sample_directories = sorted(
    [
        path
        for path in RESULT_ROOT.iterdir()
        if path.is_dir()
        and path.name.endswith(
            "_longest_isoform"
        )
    ],
    key=lambda path: path.name.casefold(),
)

if not sample_directories:
    raise SystemExit(
        "错误：没有发现*_longest_isoform结果目录"
    )

print("=" * 78)
print("复制213株BUSCO单拷贝蛋白")
print("=" * 78)
print(f"BUSCO结果：{RESULT_ROOT}")
print(f"输出目录：{OUTPUT_ROOT}")
print(
    f"发现物种结果目录："
    f"{len(sample_directories)}"
)
print()

report_rows = []
busco_presence = Counter()
sample_names = []
failed = []

for number, sample_dir in enumerate(
    sample_directories,
    start=1,
):
    sample = sample_dir.name.removesuffix(
        "_longest_isoform"
    )

    single_copy_dir = (
        sample_dir
        / "run_mucoromycota_odb12"
        / "busco_sequences"
        / "single_copy_busco_sequences"
    )

    if not single_copy_dir.is_dir():
        alternatives = list(
            sample_dir.rglob(
                "single_copy_busco_sequences"
            )
        )

        alternatives = [
            path
            for path in alternatives
            if path.is_dir()
            and "mucoromycota_odb12"
            in str(path)
        ]

        if len(alternatives) == 1:
            single_copy_dir = alternatives[0]
        else:
            failed.append(
                (
                    sample,
                    "single_copy目录不存在或不唯一",
                )
            )

            report_rows.append({
                "sample": sample,
                "source_directory": "",
                "single_copy_faa_count": 0,
                "copied_files": 0,
                "multi_record_files": 0,
                "empty_files": 0,
                "status": "MISSING_DIRECTORY",
            })

            print(
                f"[{number:03d}/"
                f"{len(sample_directories):03d}] "
                f"MISSING {sample}"
            )

            continue

    faa_files = sorted(
        single_copy_dir.glob("*.faa"),
        key=lambda path: path.name.casefold(),
    )

    species_output = (
        BY_SPECIES / sample
    )

    species_output.mkdir(
        parents=True,
        exist_ok=True,
    )

    copied = 0
    multi_record = 0
    empty_files = 0

    try:
        for source_faa in faa_files:
            busco_id = source_faa.stem

            records = read_fasta(
                source_faa
            )

            if not records:
                empty_files += 1
                continue

            if len(records) > 1:
                multi_record += 1

            species_file = (
                species_output
                / f"{busco_id}_{sample}.faa"
            )

            busco_output = (
                BY_BUSCO / busco_id
            )

            busco_output.mkdir(
                parents=True,
                exist_ok=True,
            )

            busco_species_file = (
                busco_output
                / f"{sample}.faa"
            )

            write_species_fasta(
                species_file,
                sample,
                busco_id,
                records,
            )

            shutil.copy2(
                species_file,
                busco_species_file,
            )

            copied += 1
            busco_presence[
                busco_id
            ] += 1

        sample_names.append(sample)

        status = (
            "OK"
            if copied > 0
            else "NO_FAA_FILES"
        )

        report_rows.append({
            "sample": sample,
            "source_directory": str(
                single_copy_dir
            ),
            "single_copy_faa_count": (
                len(faa_files)
            ),
            "copied_files": copied,
            "multi_record_files": (
                multi_record
            ),
            "empty_files": empty_files,
            "status": status,
        })

        print(
            f"[{number:03d}/"
            f"{len(sample_directories):03d}] "
            f"{status} {sample}: "
            f"{copied} single-copy BUSCOs"
        )

    except Exception as error:
        failed.append(
            (
                sample,
                str(error),
            )
        )

        report_rows.append({
            "sample": sample,
            "source_directory": str(
                single_copy_dir
            ),
            "single_copy_faa_count": (
                len(faa_files)
            ),
            "copied_files": copied,
            "multi_record_files": (
                multi_record
            ),
            "empty_files": empty_files,
            "status": (
                f"FAILED: {error}"
            ),
        })

        print(
            f"[{number:03d}/"
            f"{len(sample_directories):03d}] "
            f"FAILED {sample}: {error}",
            file=sys.stderr,
        )

report_columns = [
    "sample",
    "source_directory",
    "single_copy_faa_count",
    "copied_files",
    "multi_record_files",
    "empty_files",
    "status",
]

with REPORT.open(
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

with PRESENCE.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
    )

    writer.writerow([
        "BUSCO_ID",
        "species_count",
        "present_in_all_species",
    ])

    species_total = len(sample_names)

    for busco_id, count in sorted(
        busco_presence.items()
    ):
        writer.writerow([
            busco_id,
            count,
            (
                "YES"
                if count == species_total
                else "NO"
            ),
        ])

SAMPLE_LIST.write_text(
    "".join(
        f"{sample}\n"
        for sample in sorted(
            sample_names,
            key=str.casefold,
        )
    ),
    encoding="utf-8",
)

universal_count = sum(
    count == len(sample_names)
    for count in busco_presence.values()
)

print()
print("=" * 78)
print("复制完成")
print("=" * 78)
print(
    f"发现BUSCO物种目录："
    f"{len(sample_directories)}"
)
print(f"成功物种：{len(sample_names)}")
print(f"失败物种：{len(failed)}")
print(
    f"不同BUSCO marker数："
    f"{len(busco_presence)}"
)
print(
    f"所有成功物种共有marker数："
    f"{universal_count}"
)
print(f"按物种目录：{BY_SPECIES}")
print(f"按BUSCO目录：{BY_BUSCO}")
print(f"提取报告：{REPORT}")
print(f"存在统计：{PRESENCE}")
print(f"物种名单：{SAMPLE_LIST}")

if len(sample_directories) != 213:
    print(
        f"\n警告：预期213个物种目录，"
        f"实际发现{len(sample_directories)}个"
    )

if failed:
    print("\n失败物种：")

    for sample, error in failed:
        print(f"  {sample}: {error}")

    raise SystemExit(1)
