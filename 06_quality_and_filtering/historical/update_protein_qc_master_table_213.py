#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv
import re
import statistics

PROTEIN_DIR = Path(
    "/home/py/Desktop/Mucoromycota_BRAKER_core_proteomes"
)

RERUN_PROTEIN_DIR = Path(
    "/home/py/Desktop/3_rerun_protein"
)

BUSCO_MU = Path(
    "/home/py/bioinfo/mucoromycota_busco_protein_local/"
    "results_mucoromycota"
)

BUSCO_FU = Path(
    "/home/py/bioinfo/mucoromycota_busco_protein_local/"
    "results_fungi"
)

RERUN_BUSCO_MU = Path(
    "/home/py/bioinfo/mucoromycota_busco_protein_rerun3/"
    "results_mucoromycota"
)

RERUN_BUSCO_FU = Path(
    "/home/py/bioinfo/mucoromycota_busco_protein_rerun3/"
    "results_fungi"
)

BRAKER2_DIR = Path(
    "/home/py/Desktop/braker2_127_proteins"
)

BRAKER3_DIR = Path(
    "/home/py/Desktop/braker3_83_proteins"
)

OUT_CSV = Path(
    "/home/py/Desktop/protein_qc_master_table.csv"
)

OUT_TSV = Path(
    "/home/py/Desktop/protein_qc_master_table.tsv"
)

BUSCO_PATTERN = re.compile(
    r"C:\s*([0-9.]+)%\s*"
    r"\[\s*S:\s*([0-9.]+)%\s*,\s*D:\s*([0-9.]+)%\s*\]"
    r"\s*,?\s*F:\s*([0-9.]+)%"
    r"\s*,?\s*M:\s*([0-9.]+)%"
    r"\s*,?\s*n:\s*(\d+)"
)


def read_fasta_metrics(path):
    lengths = []
    ids = []
    terminal_stop = 0
    internal_stop = 0
    with_x = 0
    empty = 0

    current_id = None
    current_seq = []

    def save_sequence(seq_id, parts):
        nonlocal terminal_stop, internal_stop, with_x, empty

        if seq_id is None:
            return

        sequence = "".join(parts).replace(" ", "").upper()

        if sequence.endswith("*"):
            terminal_stop += 1

        if "*" in sequence[:-1]:
            internal_stop += 1

        if "X" in sequence:
            with_x += 1

        clean = sequence.replace("*", "").replace("-", "")
        length = len(clean)

        if length == 0:
            empty += 1

        ids.append(seq_id)
        lengths.append(length)

    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                save_sequence(current_id, current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    save_sequence(current_id, current_seq)

    valid_lengths = [x for x in lengths if x > 0]
    total_aa = sum(valid_lengths)

    if valid_lengths:
        sorted_lengths = sorted(valid_lengths, reverse=True)
        cumulative = 0
        protein_n50 = 0

        for length in sorted_lengths:
            cumulative += length
            if cumulative >= total_aa / 2:
                protein_n50 = length
                break

        mean_length = round(statistics.mean(valid_lengths), 2)
        median_length = round(statistics.median(valid_lengths), 2)
        minimum = min(valid_lengths)
        maximum = max(valid_lengths)
    else:
        protein_n50 = 0
        mean_length = 0
        median_length = 0
        minimum = 0
        maximum = 0

    duplicate_ids = len(ids) - len(set(ids))

    return {
        "protein_sequence_count": len(lengths),
        "total_aa": total_aa,
        "mean_protein_length_aa": mean_length,
        "median_protein_length_aa": median_length,
        "protein_N50_aa": protein_n50,
        "shortest_protein_aa": minimum,
        "longest_protein_aa": maximum,
        "proteins_lt_100aa": sum(x < 100 for x in valid_lengths),
        "proteins_ge_1000aa": sum(x >= 1000 for x in valid_lengths),
        "proteins_with_terminal_stop": terminal_stop,
        "proteins_with_internal_stop": internal_stop,
        "proteins_with_X": with_x,
        "empty_protein_sequences": empty,
        "duplicate_fasta_IDs": duplicate_ids,
        "protein_fasta_size_mb": round(path.stat().st_size / 1_000_000, 3),
    }


def parse_busco(busco_root, sample, lineage):
    sample_dir = busco_root / sample

    empty_result = {
        "C": "",
        "S": "",
        "D": "",
        "F": "",
        "M": "",
        "n": "",
        "single_copy_faa_count": 0,
        "status": "MISSING_RESULT_DIR",
    }

    if not sample_dir.is_dir():
        return empty_result

    summary_files = sorted(sample_dir.glob("short_summary*.json"))
    summary_files += sorted(sample_dir.glob("short_summary*.txt"))

    if not summary_files:
        summary_files = sorted(
            sample_dir.rglob("short_summary*.json")
        )
        summary_files += sorted(
            sample_dir.rglob("short_summary*.txt")
        )

    result = None

    for summary_file in summary_files:
        try:
            text = summary_file.read_text(errors="replace")
        except OSError:
            continue

        match = BUSCO_PATTERN.search(text)

        if match:
            result = {
                "C": float(match.group(1)),
                "S": float(match.group(2)),
                "D": float(match.group(3)),
                "F": float(match.group(4)),
                "M": float(match.group(5)),
                "n": int(match.group(6)),
                "status": "OK",
            }
            break

    sco_dir = (
        sample_dir
        / f"run_{lineage}"
        / "busco_sequences"
        / "single_copy_busco_sequences"
    )

    if sco_dir.is_dir():
        sco_count = sum(1 for x in sco_dir.glob("*.faa"))
    else:
        candidates = list(
            sample_dir.rglob("single_copy_busco_sequences")
        )

        if candidates:
            sco_count = sum(
                1 for x in candidates[0].glob("*.faa")
            )
        else:
            sco_count = 0

    if result is None:
        result = {
            "C": "",
            "S": "",
            "D": "",
            "F": "",
            "M": "",
            "n": "",
            "status": "SUMMARY_PARSE_FAILED",
        }

    result["single_copy_faa_count"] = sco_count
    return result


def determine_pipeline(filename):
    in_b2 = (BRAKER2_DIR / filename).is_file()
    in_b3 = (BRAKER3_DIR / filename).is_file()

    if in_b2 and in_b3:
        return "BRAKER2_AND_BRAKER3"
    if in_b2:
        return "BRAKER2"
    if in_b3:
        return "BRAKER3"

    return "UNKNOWN"


old_files = list(PROTEIN_DIR.glob("*_protein.faa"))
new_files = list(RERUN_PROTEIN_DIR.glob("*_protein.faa"))

# 新的3株优先覆盖同名旧蛋白
protein_map = {x.name: x for x in old_files}
protein_map.update({x.name: x for x in new_files})
protein_files = sorted(protein_map.values(), key=lambda x: x.name)

if not protein_files:
    raise SystemExit(
        f"ERROR: no *_protein.faa files found in {PROTEIN_DIR}"
    )

rows = []

for number, protein_file in enumerate(protein_files, 1):
    sample = protein_file.name.removesuffix("_protein.faa")

    print(
        f"[{number}/{len(protein_files)}] Processing {sample}",
        flush=True,
    )

    if protein_file.parent == RERUN_PROTEIN_DIR:
        mucor_root = RERUN_BUSCO_MU
        fungi_root = RERUN_BUSCO_FU
    else:
        mucor_root = BUSCO_MU
        fungi_root = BUSCO_FU

    mucor = parse_busco(
        mucor_root,
        sample,
        "mucoromycota_odb12",
    )

    fungi = parse_busco(
        fungi_root,
        sample,
        "fungi_odb12",
    )

    protein_metrics = read_fasta_metrics(protein_file)

    row = {
        "sample": sample,
        "annotation_pipeline": determine_pipeline(protein_file.name),
        "protein_file": protein_file.name,

        "C_mucoromycota": mucor["C"],
        "S_mucoromycota": mucor["S"],
        "D_mucoromycota": mucor["D"],
        "F_mucoromycota": mucor["F"],
        "M_mucoromycota": mucor["M"],
        "n_mucoromycota": mucor["n"],
        "single_copy_faa_mucoromycota":
            mucor["single_copy_faa_count"],
        "BUSCO_status_mucoromycota": mucor["status"],

        "C_fungi": fungi["C"],
        "S_fungi": fungi["S"],
        "D_fungi": fungi["D"],
        "F_fungi": fungi["F"],
        "M_fungi": fungi["M"],
        "n_fungi": fungi["n"],
        "single_copy_faa_fungi":
            fungi["single_copy_faa_count"],
        "BUSCO_status_fungi": fungi["status"],
    }

    row.update(protein_metrics)
    rows.append(row)


fieldnames = [
    "sample",
    "annotation_pipeline",
    "protein_file",

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


with OUT_CSV.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )
    writer.writeheader()
    writer.writerows(rows)


with OUT_TSV.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)


mu_status = Counter(
    row["BUSCO_status_mucoromycota"] for row in rows
)

fu_status = Counter(
    row["BUSCO_status_fungi"] for row in rows
)

print()
print("Protein QC master table completed")
print(f"Protein files: {len(rows)}")
print(f"CSV: {OUT_CSV}")
print(f"TSV: {OUT_TSV}")
print(f"Mucoromycota BUSCO status: {dict(mu_status)}")
print(f"Fungi BUSCO status: {dict(fu_status)}")
