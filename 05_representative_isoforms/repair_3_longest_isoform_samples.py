#!/usr/bin/env python3

from pathlib import Path
from collections import OrderedDict, defaultdict
import csv
import re
import argparse

# [release refactor] Paths are command-line arguments; repair logic unchanged.
# Historical defaults (author workstation):
#   GTF_DIR     /home/py/Desktop/Mucoromycota_BRAKER_core_gtf
#   PROTEIN_DIR /home/py/Desktop/Mucoromycota_BRAKER_core_proteomes
#   OUT_ROOT    /home/py/Desktop/Mucoromycota_BRAKER_longest_isoform
# NOTE: this repair writes INTO the OUT_ROOT produced by
# filter_longest_isoforms_from_gtf.py and replaces the files of the three
# samples listed in FAILED_SAMPLES. Its selection rule differs from the main
# script (longest CDS whose protein record exists; otherwise next candidate).
_ap = argparse.ArgumentParser(
    description="Repair longest-isoform outputs for the three historically failed samples.")
_ap.add_argument("--gtf-dir", required=True)
_ap.add_argument("--protein-dir", required=True)
_ap.add_argument("--out-root", required=True,
                 help="existing output root of filter_longest_isoforms_from_gtf.py")
_args = _ap.parse_args()

GTF_DIR = Path(_args.gtf_dir)

PROTEIN_DIR = Path(_args.protein_dir)

OUT_ROOT = Path(_args.out_root)
if not OUT_ROOT.is_dir():
    raise SystemExit(f"ERROR: --out-root does not exist: {OUT_ROOT}")

OUT_GTF = OUT_ROOT / "01_longest_isoform_gtf"
OUT_PROTEIN = OUT_ROOT / "02_longest_isoform_proteins"
OUT_MAPPING = OUT_ROOT / "03_mapping_tables"
OUT_REPORT = OUT_ROOT / "04_reports"

REPAIR_REPORT = (
    OUT_REPORT /
    "three_samples_longest_isoform_repair.tsv"
)

FAILED_SAMPLES = [
    "Linnemannia_amoeboidea",
    "Rhizopus_microsporus",
    "Rhizopus_stolonifer",
]


def parse_attributes(text):
    attributes = {}

    for key, quoted, unquoted in re.findall(
        r'([A-Za-z0-9_.:-]+)\s+'
        r'(?:"([^"]*)"|([^;\s]+))\s*;?',
        text,
    ):
        attributes[key] = quoted or unquoted

    return attributes


def line_ids(fields):
    feature = fields[2]
    attrs = parse_attributes(fields[8])

    gene_id = attrs.get("gene_id")
    transcript_id = attrs.get("transcript_id")

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
                    f"第{line_number}行不是9列GTF"
                )

            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])

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
                        f"{transcript_id}对应多个gene"
                    )

                transcript_gene[
                    transcript_id
                ] = gene_id

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
    }


def read_fasta(path):
    records = []
    header = None
    sequence = []

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
                            "".join(sequence),
                        )
                    )

                header = line[1:].strip()
                sequence = []
            else:
                if header is None:
                    raise ValueError(
                        "FASTA首个记录不是标题"
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


def fasta_candidates(header):
    candidates = [
        header.split()[0],
    ]

    for pattern in (
        r'transcript_id[=\s"]+([^";\s]+)',
        r'transcript[=\s"]+([^";\s]+)',
        r'ID[=\s"]+([^";\s]+)',
    ):
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


def choose_available_transcripts(
    data,
    protein_index,
):
    by_gene = defaultdict(list)

    for transcript_id, gene_id in (
        data["transcript_gene"].items()
    ):
        by_gene[gene_id].append(
            transcript_id
        )

    selected = {}
    audit_rows = []
    genes_without_protein = []

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
            genes_without_protein.append(
                gene_id
            )

            audit_rows.append({
                "gene_id": gene_id,
                "longest_gtf_transcript": (
                    ranked[0]
                ),
                "selected_transcript": "",
                "selection_rank": "",
                "longest_cds_bp": (
                    data["cds_length"].get(
                        ranked[0],
                        0,
                    )
                ),
                "selected_cds_bp": "",
                "number_of_isoforms": len(
                    ranked
                ),
                "status": (
                    "NO_TRANSCRIPT_HAS_PROTEIN"
                ),
            })

            continue

        chosen = available[0]
        selected[gene_id] = chosen
        selection_rank = (
            ranked.index(chosen) + 1
        )

        audit_rows.append({
            "gene_id": gene_id,
            "longest_gtf_transcript": (
                ranked[0]
            ),
            "selected_transcript": chosen,
            "selection_rank": selection_rank,
            "longest_cds_bp": (
                data["cds_length"].get(
                    ranked[0],
                    0,
                )
            ),
            "selected_cds_bp": (
                data["cds_length"].get(
                    chosen,
                    0,
                )
            ),
            "number_of_isoforms": len(
                ranked
            ),
            "status": (
                "SELECTED_LONGEST_AVAILABLE"
                if selection_rank == 1
                else "FALLBACK_TO_NEXT_AVAILABLE"
            ),
        })

    return (
        selected,
        audit_rows,
        genes_without_protein,
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

    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Derived representative-transcript GTF\n"
        )
        handle.write(
            "# For each gene, the transcript with the "
            "longest CDS that had a corresponding "
            "protein FASTA record was retained.\n"
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

            if fields[2] == "gene":
                keep = gene_id in selected_genes
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
                handle.write(line + "\n")


all_repair_rows = []

for sample in FAILED_SAMPLES:
    gtf_file = GTF_DIR / f"{sample}.gtf"
    protein_file = (
        PROTEIN_DIR /
        f"{sample}_protein.faa"
    )

    for path in (gtf_file, protein_file):
        if not path.is_file():
            raise SystemExit(
                f"文件不存在：{path}"
            )

    data = read_gtf(gtf_file)
    protein_records = read_fasta(
        protein_file
    )

    protein_index = {}

    for header, sequence in protein_records:
        for candidate in fasta_candidates(
            header
        ):
            protein_index.setdefault(
                candidate,
                (
                    header,
                    sequence,
                ),
            )

    (
        selected,
        audit_rows,
        genes_without_protein,
    ) = choose_available_transcripts(
        data,
        protein_index,
    )

    if genes_without_protein:
        print(
            f"警告：{sample}有"
            f"{len(genes_without_protein)}个gene的"
            "全部transcript都没有蛋白记录。"
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

    write_filtered_gtf(
        output_gtf,
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

    with output_protein.open(
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
        "sample",
        "gene_id",
        "longest_gtf_transcript",
        "selected_transcript",
        "selection_rank",
        "longest_cds_bp",
        "selected_cds_bp",
        "protein_id",
        "protein_length_aa",
        "number_of_isoforms",
        "status",
    ]

    audit_by_gene = {
        row["gene_id"]: row
        for row in audit_rows
    }

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

        for gene_id in sorted(
            audit_by_gene,
            key=lambda x: (
                data["gene_order"].get(
                    x,
                    10**18,
                ),
                x,
            ),
        ):
            row = audit_by_gene[gene_id]
            transcript_id = row[
                "selected_transcript"
            ]

            if transcript_id:
                header, sequence = (
                    protein_index[
                        transcript_id
                    ]
                )
                protein_id = (
                    header.split()[0]
                )
                length_aa = protein_length(
                    sequence
                )
            else:
                protein_id = ""
                length_aa = ""

            writer.writerow({
                "sample": sample,
                **row,
                "protein_id": protein_id,
                "protein_length_aa": length_aa,
            })

            all_repair_rows.append({
                "sample": sample,
                **row,
                "protein_id": protein_id,
                "protein_length_aa": length_aa,
            })

    fallbacks = sum(
        row["status"]
        == "FALLBACK_TO_NEXT_AVAILABLE"
        for row in audit_rows
    )

    print(
        f"{sample}: "
        f"代表蛋白={len(selected)}, "
        f"使用次长isoform={fallbacks}, "
        f"无可用蛋白gene="
        f"{len(genes_without_protein)}"
    )

with REPAIR_REPORT.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    columns = [
        "sample",
        "gene_id",
        "longest_gtf_transcript",
        "selected_transcript",
        "selection_rank",
        "longest_cds_bp",
        "selected_cds_bp",
        "protein_id",
        "protein_length_aa",
        "number_of_isoforms",
        "status",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(all_repair_rows)

print()
print("3株修复完成")
print(f"审计报告：{REPAIR_REPORT}")
print("原始GTF和原始蛋白均未修改。")
