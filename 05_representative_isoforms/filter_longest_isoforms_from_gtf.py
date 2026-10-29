#!/usr/bin/env python3

from pathlib import Path
from collections import OrderedDict, defaultdict
import csv
import re
import sys
import argparse

# [release refactor] Input/output locations are command-line arguments.
# Historical defaults (author workstation, not portable):
#   GTF_DIR     /home/py/Desktop/Mucoromycota_BRAKER_core_gtf
#   PROTEIN_DIR /home/py/Desktop/Mucoromycota_BRAKER_core_proteomes
#   OUT_ROOT    /home/py/Desktop/Mucoromycota_BRAKER_longest_isoform
# Selection logic below is unchanged from the historical script.
_ap = argparse.ArgumentParser(
    description="Select one representative transcript per gene from BRAKER GTF "
                "(key: total CDS length desc, transcript span desc, first GTF "
                "appearance, transcript ID) and extract matching proteins.")
_ap.add_argument("--gtf-dir", required=True, help="directory with <Strain>*.gtf from BRAKER")
_ap.add_argument("--protein-dir", required=True, help="directory with matching protein FASTA files")
_ap.add_argument("--out-root", required=True, help="output root (01_..04_ sub-directories are created)")
_ap.add_argument("--overwrite", action="store_true",
                 help="allow writing into non-empty output sub-directories (default: refuse)")
_args = _ap.parse_args()

GTF_DIR = Path(_args.gtf_dir)

PROTEIN_DIR = Path(_args.protein_dir)

OUT_ROOT = Path(_args.out_root)

OUT_GTF = OUT_ROOT / "01_longest_isoform_gtf"
OUT_PROTEIN = OUT_ROOT / "02_longest_isoform_proteins"
OUT_MAPPING = OUT_ROOT / "03_mapping_tables"
OUT_REPORT = OUT_ROOT / "04_reports"

SUMMARY = OUT_REPORT / "longest_isoform_summary.csv"
UNMATCHED = OUT_REPORT / "unmatched_transcript_ids.tsv"

for directory in (
    GTF_DIR,
    PROTEIN_DIR,
):
    if not directory.is_dir():
        raise SystemExit(
            f"错误：输入目录不存在：{directory}"
        )

if not _args.overwrite:
    for directory in (OUT_GTF, OUT_PROTEIN, OUT_MAPPING, OUT_REPORT):
        if directory.is_dir() and any(directory.iterdir()):
            raise SystemExit(
                f"错误：输出目录非空（使用 --overwrite 才允许覆盖）：{directory}"
            )

for directory in (
    OUT_GTF,
    OUT_PROTEIN,
    OUT_MAPPING,
    OUT_REPORT,
):
    directory.mkdir(parents=True, exist_ok=True)


def parse_attributes(text):
    """
    兼容常见GTF属性：
    gene_id "g1"; transcript_id "g1.t1";
    """
    attributes = {}

    for key, quoted, unquoted in re.findall(
        r'([A-Za-z0-9_.:-]+)\s+'
        r'(?:"([^"]*)"|([^;\s]+))\s*;?',
        text,
    ):
        attributes[key] = quoted or unquoted

    return attributes


def read_gtf(path):
    """
    读取GTF并计算每个transcript的CDS总长度。
    同时保存每一行，供派生GTF输出。
    """
    lines = []
    transcript_gene = {}
    cds_length = defaultdict(int)
    transcript_span = defaultdict(int)
    transcript_order = {}
    gene_order = {}
    genes_seen = set()
    transcripts_seen = set()

    with path.open(
        "r",
        encoding="utf-8",
        errors="strict",
    ) as handle:
        for line_number, raw_line in enumerate(
            handle,
            start=1,
        ):
            line = raw_line.rstrip("\n")
            lines.append(line)

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) != 9:
                raise ValueError(
                    f"第{line_number}行不是9列GTF"
                )

            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            attributes = parse_attributes(fields[8])

            gene_id = attributes.get("gene_id")
            transcript_id = attributes.get(
                "transcript_id"
            )

            # BRAKER的gene/transcript行有时第9列只有g1或g1.t1。
            if feature == "gene" and not gene_id:
                bare = fields[8].strip().strip(";")

                if (
                    bare
                    and " " not in bare
                    and '"' not in bare
                ):
                    gene_id = bare

            if feature in ("transcript", "mRNA"):
                if not transcript_id:
                    bare = fields[8].strip().strip(";")

                    if (
                        bare
                        and " " not in bare
                        and '"' not in bare
                    ):
                        transcript_id = bare

                if not gene_id and transcript_id:
                    match = re.match(
                        r"^(.*)\.t\d+$",
                        transcript_id,
                        re.I,
                    )

                    if match:
                        gene_id = match.group(1)

            if gene_id:
                genes_seen.add(gene_id)
                gene_order.setdefault(
                    gene_id,
                    line_number,
                )

            if transcript_id:
                transcripts_seen.add(transcript_id)
                transcript_order.setdefault(
                    transcript_id,
                    line_number,
                )

            if gene_id and transcript_id:
                previous = transcript_gene.get(
                    transcript_id
                )

                if (
                    previous is not None
                    and previous != gene_id
                ):
                    raise ValueError(
                        f"transcript {transcript_id} "
                        f"对应多个gene："
                        f"{previous}, {gene_id}"
                    )

                transcript_gene[transcript_id] = gene_id

            if feature == "CDS" and transcript_id:
                cds_length[transcript_id] += (
                    end - start + 1
                )

                if (
                    transcript_id
                    not in transcript_gene
                    and gene_id
                ):
                    transcript_gene[
                        transcript_id
                    ] = gene_id

            if (
                feature in ("transcript", "mRNA")
                and transcript_id
            ):
                transcript_span[transcript_id] = (
                    end - start + 1
                )

    return {
        "lines": lines,
        "transcript_gene": transcript_gene,
        "cds_length": cds_length,
        "transcript_span": transcript_span,
        "transcript_order": transcript_order,
        "gene_order": gene_order,
        "genes_seen": genes_seen,
        "transcripts_seen": transcripts_seen,
    }


def select_longest(data):
    """
    每个gene选择CDS总长度最大的transcript。
    CDS相同则依次使用：
    1. transcript跨度较长；
    2. GTF中先出现；
    3. transcript ID字典序。
    """
    by_gene = defaultdict(list)

    for transcript_id, gene_id in (
        data["transcript_gene"].items()
    ):
        by_gene[gene_id].append(
            transcript_id
        )

    selected = {}
    mappings = []

    ordered_genes = sorted(
        by_gene,
        key=lambda gene: (
            data["gene_order"].get(
                gene,
                10**18,
            ),
            gene,
        ),
    )

    for gene_id in ordered_genes:
        transcripts = by_gene[gene_id]

        ranked = sorted(
            transcripts,
            key=lambda transcript: (
                -data["cds_length"].get(
                    transcript,
                    0,
                ),
                -data["transcript_span"].get(
                    transcript,
                    0,
                ),
                data["transcript_order"].get(
                    transcript,
                    10**18,
                ),
                transcript,
            ),
        )

        chosen = ranked[0]
        selected[gene_id] = chosen

        mappings.append({
            "gene_id": gene_id,
            "selected_transcript_id": chosen,
            "selected_cds_length_bp": (
                data["cds_length"].get(
                    chosen,
                    0,
                )
            ),
            "selected_transcript_span_bp": (
                data["transcript_span"].get(
                    chosen,
                    0,
                )
            ),
            "number_of_isoforms": len(
                transcripts
            ),
            "all_transcript_ids": ",".join(
                sorted(transcripts)
            ),
        })

    return selected, mappings


def read_fasta(path):
    records = []
    header = None
    sequence_parts = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="strict",
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
                        "FASTA第一个有效记录不是标题"
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


def fasta_candidates(header):
    """
    从BRAKER蛋白标题中提取可能的transcript ID。
    第一个token通常就是g1.t1。
    也兼容 transcript_id= 和 transcript_id 字段。
    """
    candidates = []
    first = header.split()[0]

    candidates.extend([
        first,
        first.lstrip(">"),
    ])

    patterns = [
        r'transcript_id[=\s"]+([^";\s]+)',
        r'transcript[=\s"]+([^";\s]+)',
        r'ID[=\s"]+([^";\s]+)',
    ]

    for pattern in patterns:
        for match in re.findall(
            pattern,
            header,
            re.I,
        ):
            candidates.append(match)

    # BRAKER偶尔在标题中加入gene前缀或蛋白后缀。
    expanded = []

    for candidate in candidates:
        candidate = candidate.strip(
            '";|,'
        )

        expanded.append(candidate)

        for suffix in (
            ".aa",
            ".pep",
            ".protein",
        ):
            if candidate.lower().endswith(
                suffix
            ):
                expanded.append(
                    candidate[:-len(suffix)]
                )

    return list(
        OrderedDict.fromkeys(expanded)
    )


def protein_length(sequence):
    return len(
        sequence.rstrip("*").replace("-", "")
    )


def wrap_sequence(sequence, width=60):
    return "\n".join(
        sequence[index:index + width]
        for index in range(
            0,
            len(sequence),
            width,
        )
    )


def line_ids(fields):
    """
    返回当前GTF行的gene ID和transcript ID。
    """
    feature = fields[2]
    attrs = parse_attributes(fields[8])
    gene_id = attrs.get("gene_id")
    transcript_id = attrs.get(
        "transcript_id"
    )

    if feature == "gene" and not gene_id:
        bare = fields[8].strip().strip(";")

        if (
            bare
            and " " not in bare
            and '"' not in bare
        ):
            gene_id = bare

    if feature in ("transcript", "mRNA"):
        if not transcript_id:
            bare = fields[8].strip().strip(";")

            if (
                bare
                and " " not in bare
                and '"' not in bare
            ):
                transcript_id = bare

        if not gene_id and transcript_id:
            match = re.match(
                r"^(.*)\.t\d+$",
                transcript_id,
                re.I,
            )

            if match:
                gene_id = match.group(1)

    return gene_id, transcript_id


def write_filtered_gtf(
    output_path,
    data,
    selected,
):
    """
    生成派生的最长isoform GTF。
    原始GTF只读，不执行任何写入。
    """
    selected_transcripts = set(
        selected.values()
    )
    selected_genes = set(
        selected.keys()
    )

    kept = 0

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Derived representative-transcript GTF\n"
        )
        handle.write(
            "# One transcript with the longest total "
            "CDS length was retained per gene.\n"
        )

        for line in data["lines"]:
            if not line:
                continue

            if line.startswith("#"):
                handle.write(line + "\n")
                continue

            fields = line.split("\t")
            gene_id, transcript_id = line_ids(
                fields
            )

            keep = False

            if fields[2] == "gene":
                keep = (
                    gene_id in selected_genes
                )

            elif transcript_id:
                keep = (
                    transcript_id
                    in selected_transcripts
                )

            elif gene_id:
                keep = (
                    gene_id in selected_genes
                )

            if keep:
                handle.write(line + "\n")
                kept += 1

    return kept


gtf_files = sorted(
    GTF_DIR.glob("*.gtf"),
    key=lambda path: path.name.casefold(),
)

protein_files = {
    path.name.removesuffix(
        "_protein.faa"
    ): path
    for path in PROTEIN_DIR.glob(
        "*_protein.faa"
    )
}

if not gtf_files:
    raise SystemExit("错误：没有发现GTF文件")

print("=" * 78)
print("从GTF选择每个gene的最长编码isoform")
print("=" * 78)
print(f"GTF文件数：{len(gtf_files)}")
print(f"蛋白文件数：{len(protein_files)}")
print("原始GTF：只读，绝不修改")
print()

summary_rows = []
unmatched_rows = []
failed = []

for number, gtf_file in enumerate(
    gtf_files,
    start=1,
):
    sample = gtf_file.stem
    protein_file = protein_files.get(sample)

    status = "OK"

    try:
        if protein_file is None:
            raise ValueError(
                "没有同名的 "
                f"{sample}_protein.faa"
            )

        data = read_gtf(gtf_file)
        selected, mappings = select_longest(
            data
        )

        if not selected:
            raise ValueError(
                "没有从GTF解析到gene-transcript关系"
            )

        protein_records = read_fasta(
            protein_file
        )

        protein_index = {}
        duplicate_fasta_ids = set()

        for header, sequence in protein_records:
            matched_ids = fasta_candidates(
                header
            )

            for candidate in matched_ids:
                if candidate in protein_index:
                    duplicate_fasta_ids.add(
                        candidate
                    )
                else:
                    protein_index[candidate] = (
                        header,
                        sequence,
                    )

        selected_records = []
        missing_ids = []
        selected_map = {}

        for gene_id, transcript_id in (
            selected.items()
        ):
            record = protein_index.get(
                transcript_id
            )

            if record is None:
                missing_ids.append(
                    transcript_id
                )

                unmatched_rows.append({
                    "sample": sample,
                    "gene_id": gene_id,
                    "transcript_id": (
                        transcript_id
                    ),
                    "problem": (
                        "selected_transcript_not_"
                        "found_in_protein_fasta"
                    ),
                })

                continue

            header, sequence = record

            selected_records.append(
                (
                    gene_id,
                    transcript_id,
                    header,
                    sequence,
                )
            )

            selected_map[transcript_id] = {
                "protein_id": (
                    header.split()[0]
                ),
                "protein_length_aa": (
                    protein_length(sequence)
                ),
            }

        # 任何缺失ID都不输出半成品，避免误用于BUSCO。
        if missing_ids:
            raise ValueError(
                f"{len(missing_ids)}个所选transcript "
                "在蛋白FASTA中找不到"
            )

        output_gtf = (
            OUT_GTF /
            f"{sample}_longest_isoform.gtf"
        )

        output_protein = (
            OUT_PROTEIN /
            f"{sample}_longest_isoform.faa"
        )

        mapping_file = (
            OUT_MAPPING /
            f"{sample}_longest_isoform_mapping.tsv"
        )

        kept_gtf_lines = write_filtered_gtf(
            output_gtf,
            data,
            selected,
        )

        selected_records.sort(
            key=lambda record: (
                data["gene_order"].get(
                    record[0],
                    10**18,
                ),
                record[0],
            )
        )

        with output_protein.open(
            "w",
            encoding="utf-8",
        ) as handle:
            for (
                gene_id,
                transcript_id,
                header,
                sequence,
            ) in selected_records:
                handle.write(
                    f">{header}\n"
                )
                handle.write(
                    wrap_sequence(sequence)
                    + "\n"
                )

        mapping_columns = [
            "sample",
            "gene_id",
            "selected_transcript_id",
            "selected_protein_id",
            "selected_cds_length_bp",
            "protein_length_aa",
            "selected_transcript_span_bp",
            "number_of_isoforms",
            "all_transcript_ids",
        ]

        with mapping_file.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=mapping_columns,
                delimiter="\t",
            )
            writer.writeheader()

            for mapping in mappings:
                transcript_id = mapping[
                    "selected_transcript_id"
                ]

                protein_info = selected_map[
                    transcript_id
                ]

                writer.writerow({
                    "sample": sample,
                    "gene_id": (
                        mapping["gene_id"]
                    ),
                    "selected_transcript_id": (
                        transcript_id
                    ),
                    "selected_protein_id": (
                        protein_info["protein_id"]
                    ),
                    "selected_cds_length_bp": (
                        mapping[
                            "selected_cds_length_bp"
                        ]
                    ),
                    "protein_length_aa": (
                        protein_info[
                            "protein_length_aa"
                        ]
                    ),
                    "selected_transcript_span_bp": (
                        mapping[
                            "selected_transcript_span_bp"
                        ]
                    ),
                    "number_of_isoforms": (
                        mapping[
                            "number_of_isoforms"
                        ]
                    ),
                    "all_transcript_ids": (
                        mapping[
                            "all_transcript_ids"
                        ]
                    ),
                })

        genes_with_multiple = sum(
            int(
                mapping["number_of_isoforms"]
            ) > 1
            for mapping in mappings
        )

        summary_rows.append({
            "sample": sample,
            "original_gtf": str(gtf_file),
            "original_protein": str(
                protein_file
            ),
            "gtf_genes": len(selected),
            "gtf_transcripts": len(
                data["transcript_gene"]
            ),
            "original_proteins": len(
                protein_records
            ),
            "representative_proteins": len(
                selected_records
            ),
            "removed_isoforms": (
                len(protein_records)
                - len(selected_records)
            ),
            "genes_with_multiple_isoforms": (
                genes_with_multiple
            ),
            "selected_ids_missing_in_fasta": 0,
            "duplicate_fasta_ids": len(
                duplicate_fasta_ids
            ),
            "filtered_gtf_data_lines": (
                kept_gtf_lines
            ),
            "status": status,
        })

        print(
            f"[{number:03d}/{len(gtf_files):03d}] "
            f"OK {sample}: "
            f"transcripts="
            f"{len(data['transcript_gene'])}, "
            f"representatives="
            f"{len(selected_records)}, "
            f"removed="
            f"{len(protein_records) - len(selected_records)}"
        )

    except Exception as error:
        failed.append(
            (
                sample,
                str(error),
            )
        )

        # 删除本次可能产生的派生半成品。
        for possible_output in (
            OUT_GTF /
            f"{sample}_longest_isoform.gtf",
            OUT_PROTEIN /
            f"{sample}_longest_isoform.faa",
            OUT_MAPPING /
            f"{sample}_longest_isoform_mapping.tsv",
        ):
            if possible_output.exists():
                possible_output.unlink()

        summary_rows.append({
            "sample": sample,
            "original_gtf": str(gtf_file),
            "original_protein": (
                str(protein_file)
                if protein_file
                else ""
            ),
            "gtf_genes": "",
            "gtf_transcripts": "",
            "original_proteins": "",
            "representative_proteins": "",
            "removed_isoforms": "",
            "genes_with_multiple_isoforms": "",
            "selected_ids_missing_in_fasta": "",
            "duplicate_fasta_ids": "",
            "filtered_gtf_data_lines": "",
            "status": f"FAILED: {error}",
        })

        print(
            f"[{number:03d}/{len(gtf_files):03d}] "
            f"FAILED {sample}: {error}",
            file=sys.stderr,
        )

summary_columns = [
    "sample",
    "original_gtf",
    "original_protein",
    "gtf_genes",
    "gtf_transcripts",
    "original_proteins",
    "representative_proteins",
    "removed_isoforms",
    "genes_with_multiple_isoforms",
    "selected_ids_missing_in_fasta",
    "duplicate_fasta_ids",
    "filtered_gtf_data_lines",
    "status",
]

with SUMMARY.open(
    "w",
    encoding="utf-8-sig",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_columns,
    )
    writer.writeheader()
    writer.writerows(summary_rows)

unmatched_columns = [
    "sample",
    "gene_id",
    "transcript_id",
    "problem",
]

with UNMATCHED.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=unmatched_columns,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(unmatched_rows)

print()
print("=" * 78)
print("处理完成")
print("=" * 78)
print(f"GTF输入：{len(gtf_files)}")
print(
    f"成功："
    f"{len(gtf_files) - len(failed)}"
)
print(f"失败：{len(failed)}")
print(f"输出根目录：{OUT_ROOT}")
print(f"汇总报告：{SUMMARY}")
print(f"未匹配ID：{UNMATCHED}")
print()
print("原始GTF和原始蛋白文件均未修改。")

if failed:
    print("\n失败样本：")

    for sample, error in failed:
        print(f"  {sample}: {error}")

    raise SystemExit(1)
