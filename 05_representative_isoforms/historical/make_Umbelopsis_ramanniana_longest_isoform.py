#!/usr/bin/env python3

from pathlib import Path
from collections import OrderedDict, defaultdict
import csv
import re
import shutil

WORK_DIR = Path(
    "/home/py/Desktop/Umbelopsis_ramanniana"
)

GTF_FILE = (
    WORK_DIR /
    "Umbelopsis_ramanniana.gtf"
)

PROTEIN_FILE = (
    WORK_DIR /
    "Umbelopsis_ramanniana_protein.faa"
)

OUTPUT_GTF = (
    WORK_DIR /
    "Umbelopsis_ramanniana_longest_isoform.gtf"
)

OUTPUT_PROTEIN = (
    WORK_DIR /
    "Umbelopsis_ramanniana_longest_isoform.faa"
)

MAPPING_FILE = (
    WORK_DIR /
    "Umbelopsis_ramanniana_longest_isoform_mapping.tsv"
)

UNMATCHED_FILE = (
    WORK_DIR /
    "Umbelopsis_ramanniana_unmatched_genes.tsv"
)

SUMMARY_FILE = (
    WORK_DIR /
    "Umbelopsis_ramanniana_longest_isoform_summary.txt"
)

for path in (GTF_FILE, PROTEIN_FILE):
    if not path.is_file():
        raise SystemExit(
            f"错误：输入文件不存在：{path}"
        )


def parse_attributes(text):
    attributes = {}

    for key, quoted, unquoted in re.findall(
        r'([A-Za-z0-9_.:-]+)\s+'
        r'(?:"([^"]*)"|([^;\s]+))\s*;?',
        text,
    ):
        attributes[key] = (
            quoted or unquoted
        )

    return attributes


def line_ids(fields):
    feature = fields[2]
    attributes = parse_attributes(
        fields[8]
    )

    gene_id = attributes.get("gene_id")
    transcript_id = attributes.get(
        "transcript_id"
    )

    # 兼容BRAKER的gene行：
    # 第9列可能只有g1，而不是gene_id "g1";
    if feature == "gene" and not gene_id:
        bare = fields[8].strip().strip(";")

        if (
            bare
            and " " not in bare
            and '"' not in bare
        ):
            gene_id = bare

    # 兼容BRAKER的transcript行：
    # 第9列可能只有g1.t1。
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


def read_gtf(path):
    lines = []
    transcript_gene = {}
    cds_length = defaultdict(int)
    transcript_span = defaultdict(int)
    transcript_order = {}
    gene_order = {}

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
                    f"GTF第{line_number}行不是9列"
                )

            feature = fields[2]

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                raise ValueError(
                    f"GTF第{line_number}行坐标不是整数"
                )

            gene_id, transcript_id = line_ids(
                fields
            )

            if gene_id:
                gene_order.setdefault(
                    gene_id,
                    line_number,
                )

            if transcript_id:
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

                transcript_gene[
                    transcript_id
                ] = gene_id

            if feature == "CDS" and transcript_id:
                cds_length[
                    transcript_id
                ] += end - start + 1

                if (
                    transcript_id
                    not in transcript_gene
                    and gene_id
                ):
                    transcript_gene[
                        transcript_id
                    ] = gene_id

            if (
                feature in (
                    "transcript",
                    "mRNA",
                )
                and transcript_id
            ):
                transcript_span[
                    transcript_id
                ] = end - start + 1

    return {
        "lines": lines,
        "transcript_gene": transcript_gene,
        "cds_length": cds_length,
        "transcript_span": transcript_span,
        "transcript_order": transcript_order,
        "gene_order": gene_order,
    }


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
    candidates = [
        header.split()[0],
    ]

    patterns = [
        r'transcript_id[=\s"]+([^";\s]+)',
        r'transcript[=\s"]+([^";\s]+)',
        r'ID[=\s"]+([^";\s]+)',
    ]

    for pattern in patterns:
        candidates.extend(
            re.findall(
                pattern,
                header,
                re.I,
            )
        )

    cleaned = []

    for candidate in candidates:
        candidate = candidate.strip(
            '";|,'
        )

        cleaned.append(candidate)

        for suffix in (
            ".aa",
            ".pep",
            ".protein",
        ):
            if candidate.lower().endswith(
                suffix
            ):
                cleaned.append(
                    candidate[:-len(suffix)]
                )

    return list(
        OrderedDict.fromkeys(cleaned)
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


def choose_longest_available(
    data,
    protein_index,
):
    """
    每个gene按以下顺序选择：
    1. CDS总长度最长；
    2. transcript跨度最长；
    3. GTF中先出现；
    4. transcript ID字典序。

    如果最长转录本不在蛋白FASTA中，
    则选择下一个有蛋白记录的最长转录本。
    """
    by_gene = defaultdict(list)

    for transcript_id, gene_id in (
        data["transcript_gene"].items()
    ):
        by_gene[gene_id].append(
            transcript_id
        )

    selected = {}
    mapping_rows = []
    unmatched_rows = []

    ordered_genes = sorted(
        by_gene,
        key=lambda gene_id: (
            data["gene_order"].get(
                gene_id,
                10**18,
            ),
            gene_id,
        ),
    )

    for gene_id in ordered_genes:
        ranked = sorted(
            by_gene[gene_id],
            key=lambda transcript_id: (
                -data["cds_length"].get(
                    transcript_id,
                    0,
                ),
                -data["transcript_span"].get(
                    transcript_id,
                    0,
                ),
                data["transcript_order"].get(
                    transcript_id,
                    10**18,
                ),
                transcript_id,
            ),
        )

        available = [
            transcript_id
            for transcript_id in ranked
            if transcript_id in protein_index
        ]

        if not available:
            unmatched_rows.append({
                "gene_id": gene_id,
                "longest_gtf_transcript": (
                    ranked[0]
                ),
                "number_of_isoforms": (
                    len(ranked)
                ),
                "all_transcript_ids": (
                    ",".join(ranked)
                ),
                "reason": (
                    "NO_TRANSCRIPT_HAS_"
                    "PROTEIN_RECORD"
                ),
            })
            continue

        chosen = available[0]
        selection_rank = (
            ranked.index(chosen) + 1
        )

        selected[gene_id] = chosen

        header, sequence = protein_index[
            chosen
        ]

        mapping_rows.append({
            "gene_id": gene_id,
            "longest_gtf_transcript": (
                ranked[0]
            ),
            "selected_transcript_id": (
                chosen
            ),
            "selection_rank": (
                selection_rank
            ),
            "selected_protein_id": (
                header.split()[0]
            ),
            "selected_cds_length_bp": (
                data["cds_length"].get(
                    chosen,
                    0,
                )
            ),
            "selected_protein_length_aa": (
                protein_length(sequence)
            ),
            "number_of_isoforms": (
                len(ranked)
            ),
            "all_transcript_ids": (
                ",".join(ranked)
            ),
            "status": (
                "SELECTED_LONGEST"
                if selection_rank == 1
                else
                "FALLBACK_TO_NEXT_AVAILABLE"
            ),
        })

    return (
        selected,
        mapping_rows,
        unmatched_rows,
    )


def write_filtered_gtf(
    path,
    data,
    selected,
):
    selected_genes = set(selected)
    selected_transcripts = set(
        selected.values()
    )

    kept_data_lines = 0

    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Derived representative-transcript GTF\n"
        )
        handle.write(
            "# One transcript with the longest "
            "total CDS length and an available "
            "protein record was retained per gene.\n"
        )

        for line in data["lines"]:
            if not line:
                continue

            if line.startswith("#"):
                handle.write(
                    line + "\n"
                )
                continue

            fields = line.split("\t")
            gene_id, transcript_id = line_ids(
                fields
            )

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
            else:
                keep = False

            if keep:
                handle.write(
                    line + "\n"
                )
                kept_data_lines += 1

    return kept_data_lines


print("读取GTF……")
data = read_gtf(GTF_FILE)

if not data["transcript_gene"]:
    raise SystemExit(
        "错误：没有从GTF解析到gene-transcript关系"
    )

print("读取蛋白FASTA……")
protein_records = read_fasta(
    PROTEIN_FILE
)

if not protein_records:
    raise SystemExit(
        "错误：蛋白FASTA中没有记录"
    )

protein_index = {}
duplicate_index_ids = set()

for header, sequence in protein_records:
    for candidate in fasta_candidates(
        header
    ):
        if candidate in protein_index:
            duplicate_index_ids.add(
                candidate
            )
        else:
            protein_index[candidate] = (
                header,
                sequence,
            )

(
    selected,
    mapping_rows,
    unmatched_rows,
) = choose_longest_available(
    data,
    protein_index,
)

if not selected:
    raise SystemExit(
        "错误：没有任何GTF transcript与蛋白FASTA匹配"
    )

# 正式写入前备份已有派生结果。
for output in (
    OUTPUT_GTF,
    OUTPUT_PROTEIN,
    MAPPING_FILE,
    UNMATCHED_FILE,
    SUMMARY_FILE,
):
    if output.exists():
        backup = output.with_suffix(
            output.suffix + ".previous"
        )

        shutil.copy2(
            output,
            backup,
        )

kept_gtf_lines = write_filtered_gtf(
    OUTPUT_GTF,
    data,
    selected,
)

ordered_genes = sorted(
    selected,
    key=lambda gene_id: (
        data["gene_order"].get(
            gene_id,
            10**18,
        ),
        gene_id,
    ),
)

with OUTPUT_PROTEIN.open(
    "w",
    encoding="utf-8",
) as handle:
    for gene_id in ordered_genes:
        transcript_id = selected[
            gene_id
        ]

        header, sequence = protein_index[
            transcript_id
        ]

        handle.write(f">{header}\n")
        handle.write(
            wrap_sequence(sequence) + "\n"
        )

mapping_columns = [
    "gene_id",
    "longest_gtf_transcript",
    "selected_transcript_id",
    "selection_rank",
    "selected_protein_id",
    "selected_cds_length_bp",
    "selected_protein_length_aa",
    "number_of_isoforms",
    "all_transcript_ids",
    "status",
]

with MAPPING_FILE.open(
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
    writer.writerows(mapping_rows)

unmatched_columns = [
    "gene_id",
    "longest_gtf_transcript",
    "number_of_isoforms",
    "all_transcript_ids",
    "reason",
]

with UNMATCHED_FILE.open(
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

fallback_count = sum(
    row["status"]
    == "FALLBACK_TO_NEXT_AVAILABLE"
    for row in mapping_rows
)

multi_isoform_genes = sum(
    row["number_of_isoforms"] > 1
    for row in mapping_rows
)

summary = (
    "Umbelopsis_ramanniana longest-isoform summary\n"
    "==================================================\n"
    f"Original GTF: {GTF_FILE}\n"
    f"Original protein FASTA: {PROTEIN_FILE}\n"
    f"Original protein records: {len(protein_records)}\n"
    f"GTF transcripts: {len(data['transcript_gene'])}\n"
    f"Selected representative proteins: {len(selected)}\n"
    f"Removed alternative protein records: "
    f"{len(protein_records) - len(selected)}\n"
    f"Genes with multiple isoforms: "
    f"{multi_isoform_genes}\n"
    f"Fallback to next available isoform: "
    f"{fallback_count}\n"
    f"Genes without any protein record: "
    f"{len(unmatched_rows)}\n"
    f"Duplicate FASTA index IDs: "
    f"{len(duplicate_index_ids)}\n"
    f"Derived GTF data lines: "
    f"{kept_gtf_lines}\n"
    f"Output GTF: {OUTPUT_GTF}\n"
    f"Output protein: {OUTPUT_PROTEIN}\n"
)

SUMMARY_FILE.write_text(
    summary,
    encoding="utf-8",
)

print()
print(summary)

if unmatched_rows:
    print(
        "注意：存在GTF gene没有任何对应蛋白，"
        "它们保留在原始GTF中，但未进入派生结果。"
    )

print(
    "完成：原始GTF和原始蛋白均未修改。"
)
