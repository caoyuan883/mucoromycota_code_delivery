#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import bz2
import csv
import gzip
import io
import lzma
import re
import shutil
import tarfile
import zipfile
import argparse

# [release refactor] paths are command-line arguments; extraction logic unchanged.
# Historical values (author workstation):
#   SOURCE_DIRS /home/py/Desktop/braker2_127_core_results
#               /home/py/Desktop/braker_3_86_core_results/_delivery
#   OUTPUT_DIR  /home/py/Desktop/Mucoromycota_BRAKER_core_gtf
_ap = argparse.ArgumentParser(description="Extract braker.gtf from packaged BRAKER2/BRAKER3 core-result archives.")
_ap.add_argument("--source-dirs", nargs="+", required=True,
                 help="archive folders in historical order: BRAKER2 core results, then BRAKER3 _delivery")
_ap.add_argument("--out", required=True, help="output folder (must not exist or be empty)")
_args = _ap.parse_args()

SOURCE_DIRS = [Path(p) for p in _args.source_dirs]
for _d in SOURCE_DIRS:
    if not _d.is_dir():
        raise SystemExit(f"ERROR: source folder not found: {_d}")

OUTPUT_DIR = Path(_args.out)
if OUTPUT_DIR.is_dir() and any(OUTPUT_DIR.iterdir()):
    raise SystemExit(f"ERROR: output folder not empty: {OUTPUT_DIR}")

REPORT = OUTPUT_DIR / "gtf_extraction_report.tsv"

ARCHIVE_SUFFIXES = (
    ".tar.gz",
    ".tgz",
    ".tar.bz2",
    ".tbz2",
    ".tar.xz",
    ".txz",
    ".tar",
    ".zip",
    ".gz",
    ".bz2",
    ".xz",
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def archive_type(path):
    name = path.name.lower()

    if name.endswith(".zip"):
        return "zip"

    if name.endswith(
        (
            ".tar.gz",
            ".tgz",
            ".tar.bz2",
            ".tbz2",
            ".tar.xz",
            ".txz",
            ".tar",
        )
    ):
        return "tar"

    if name.endswith(".gz"):
        return "gzip"

    if name.endswith(".bz2"):
        return "bz2"

    if name.endswith(".xz"):
        return "xz"

    return None


def remove_archive_suffix(filename):
    lower = filename.lower()

    for suffix in sorted(
        ARCHIVE_SUFFIXES,
        key=len,
        reverse=True,
    ):
        if lower.endswith(suffix):
            return filename[:-len(suffix)]

    return filename


def clean_sample_name(archive):
    """
    从压缩包文件名取得菌株名，并删除常见交付后缀。
    """
    sample = remove_archive_suffix(archive.name)

    removable_suffixes = (
        "_core_results",
        "_core_result",
        "_delivery",
        "_braker_results",
        "_braker_result",
        "_braker3_results",
        "_braker3_result",
        "_braker2_results",
        "_braker2_result",
        "_results",
        "_result",
        "_braker3",
        "_braker2",
        "_braker",
        ".braker3",
        ".braker2",
        ".braker",
    )

    changed = True

    while changed:
        changed = False

        for suffix in removable_suffixes:
            if sample.lower().endswith(suffix.lower()):
                sample = sample[:-len(suffix)]
                changed = True
                break

    return sample.strip("._- ")


def basename(member_name):
    return member_name.replace("\\", "/").rstrip("/").split("/")[-1]


def member_priority(member_name):
    """
    数字越小，优先级越高。
    只接受精确文件名，防止误选Augustus或GeneMark内部GTF。
    """
    base = basename(member_name).lower()

    if base == "braker_final.gtf":
        return 1

    if base == "braker.gtf":
        return 2

    return None


def choose_member(member_names):
    candidates = []

    for member_name in member_names:
        priority = member_priority(member_name)

        if priority is not None:
            candidates.append(
                (
                    priority,
                    member_name.count("/")
                    + member_name.count("\\"),
                    len(member_name),
                    member_name,
                )
            )

    if not candidates:
        return None, []

    candidates.sort()
    selected = candidates[0][3]

    return selected, [item[3] for item in candidates]


def extract_from_tar(archive):
    with tarfile.open(archive, "r:*") as tar:
        files = [
            member
            for member in tar.getmembers()
            if member.isfile()
        ]

        selected_name, candidates = choose_member(
            [member.name for member in files]
        )

        if selected_name is None:
            return None, None, candidates

        member = next(
            item
            for item in files
            if item.name == selected_name
        )

        source = tar.extractfile(member)

        if source is None:
            raise RuntimeError(
                f"无法读取成员：{selected_name}"
            )

        content = source.read()

        return content, selected_name, candidates


def extract_from_zip(archive):
    with zipfile.ZipFile(archive, "r") as zf:
        files = [
            info
            for info in zf.infolist()
            if not info.is_dir()
        ]

        selected_name, candidates = choose_member(
            [info.filename for info in files]
        )

        if selected_name is None:
            return None, None, candidates

        content = zf.read(selected_name)

        return content, selected_name, candidates


def extract_single_compressed(archive, compression):
    """
    兼容单个braker.gtf.gz、braker_final.gtf.gz等文件。
    """
    decompressed_name = remove_archive_suffix(
        archive.name
    )

    if member_priority(decompressed_name) is None:
        return None, None, []

    if compression == "gzip":
        opener = gzip.open
    elif compression == "bz2":
        opener = bz2.open
    elif compression == "xz":
        opener = lzma.open
    else:
        raise ValueError(compression)

    with opener(archive, "rb") as handle:
        return (
            handle.read(),
            decompressed_name,
            [decompressed_name],
        )


def validate_gtf(content):
    """
    做最低限度检查，避免输出空文件或HTML错误页。
    """
    if not content:
        return False, "empty_file"

    text = content.decode(
        "utf-8",
        errors="replace",
    )

    valid_features = {
        "gene",
        "transcript",
        "mRNA",
        "exon",
        "CDS",
        "start_codon",
        "stop_codon",
    }

    valid_line_count = 0

    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")

        if (
            len(fields) >= 9
            and fields[2] in valid_features
        ):
            valid_line_count += 1

            if valid_line_count >= 1:
                return True, "OK"

    return False, "no_valid_gtf_records"


archives = []

for source_dir in SOURCE_DIRS:
    if not source_dir.is_dir():
        print(
            f"警告：源目录不存在，跳过：{source_dir}"
        )
        continue

    for path in source_dir.rglob("*"):
        if path.is_file() and archive_type(path):
            archives.append(
                (
                    source_dir,
                    path,
                )
            )

archives.sort(
    key=lambda item: (
        str(item[0]).casefold(),
        item[1].name.casefold(),
        str(item[1]).casefold(),
    )
)

if not archives:
    raise SystemExit(
        "错误：两个源目录中均未发现支持的压缩包"
    )

print("=" * 78)
print("从BRAKER压缩包复制最终GTF")
print("=" * 78)
print(f"找到压缩包：{len(archives)}")
print(f"输出目录：{OUTPUT_DIR}")
print()

rows = []
successful_outputs = {}
counts = Counter()

for number, (source_dir, archive) in enumerate(
    archives,
    start=1,
):
    sample = clean_sample_name(archive)
    kind = archive_type(archive)

    if not sample:
        sample = remove_archive_suffix(
            archive.name
        )

    output = (
        OUTPUT_DIR /
        f"{sample}_braker.gtf"
    )

    row = {
        "sample": sample,
        "source_group": source_dir.name,
        "archive": str(archive),
        "archive_type": kind or "",
        "selected_member": "",
        "candidate_count": 0,
        "output": str(output),
        "output_size_bytes": "",
        "status": "",
    }

    try:
        if kind == "tar":
            content, member, candidates = (
                extract_from_tar(archive)
            )
        elif kind == "zip":
            content, member, candidates = (
                extract_from_zip(archive)
            )
        elif kind in {"gzip", "bz2", "xz"}:
            content, member, candidates = (
                extract_single_compressed(
                    archive,
                    kind,
                )
            )
        else:
            raise RuntimeError(
                "不支持的压缩格式"
            )

        row["candidate_count"] = len(candidates)

        if content is None:
            row["status"] = (
                "MISSING: no braker_final.gtf "
                "or braker.gtf"
            )
            counts["missing"] += 1

            print(
                f"[{number:03d}/{len(archives):03d}] "
                f"MISSING  {archive.name}"
            )

            rows.append(row)
            continue

        valid, reason = validate_gtf(content)

        if not valid:
            row["selected_member"] = member or ""
            row["status"] = f"INVALID: {reason}"
            counts["invalid"] += 1

            print(
                f"[{number:03d}/{len(archives):03d}] "
                f"INVALID  {archive.name}: {reason}"
            )

            rows.append(row)
            continue

        if sample in successful_outputs:
            previous = successful_outputs[sample]

            row["selected_member"] = member or ""
            row["status"] = (
                "DUPLICATE_SAMPLE_NOT_WRITTEN: "
                f"already extracted from {previous}"
            )
            counts["duplicate_sample"] += 1

            print(
                f"[{number:03d}/{len(archives):03d}] "
                f"DUPLICATE {sample}: "
                f"未覆盖已有输出"
            )

            rows.append(row)
            continue

        # 直接把压缩包成员内容写到新文件。
        # 不对原压缩包执行写入或删除操作。
        with output.open("wb") as handle:
            handle.write(content)

        successful_outputs[sample] = archive
        row["selected_member"] = member or ""
        row["output_size_bytes"] = output.stat().st_size
        row["status"] = "OK"

        selected_base = basename(member).lower()

        if selected_base == "braker_final.gtf":
            counts["braker_final.gtf"] += 1
        elif selected_base == "braker.gtf":
            counts["braker.gtf"] += 1

        counts["ok"] += 1

        print(
            f"[{number:03d}/{len(archives):03d}] "
            f"OK  {sample} <- {member}"
        )

    except Exception as error:
        row["status"] = (
            f"FAILED: {type(error).__name__}: "
            f"{error}"
        )
        counts["failed"] += 1

        print(
            f"[{number:03d}/{len(archives):03d}] "
            f"FAILED {archive.name}: {error}"
        )

    rows.append(row)

columns = [
    "sample",
    "source_group",
    "archive",
    "archive_type",
    "selected_member",
    "candidate_count",
    "output",
    "output_size_bytes",
    "status",
]

with REPORT.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

output_gtfs = sorted(
    OUTPUT_DIR.glob("*_braker.gtf")
)

print()
print("=" * 78)
print("提取完成")
print("=" * 78)
print(f"扫描压缩包数：{len(archives)}")
print(f"成功输出GTF：{len(output_gtfs)}")
print(
    f"其中braker_final.gtf："
    f"{counts['braker_final.gtf']}"
)
print(
    f"其中braker.gtf："
    f"{counts['braker.gtf']}"
)
print(f"未找到GTF：{counts['missing']}")
print(f"无效GTF：{counts['invalid']}")
print(
    f"同名菌株未覆盖："
    f"{counts['duplicate_sample']}"
)
print(f"提取失败：{counts['failed']}")
print(f"输出目录：{OUTPUT_DIR}")
print(f"报告文件：{REPORT}")
