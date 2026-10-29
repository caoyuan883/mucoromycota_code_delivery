#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
from pathlib import Path

COLUMNS = [
    "sample", "status", "selected_clade", "closest_species", "closest_taxid",
    "conserved_hogs", "proteins_total", "complete_n", "complete_pct",
    "single_n", "single_pct", "duplicated_n", "duplicated_pct",
    "duplicated_unexpected_n", "duplicated_unexpected_pct",
    "duplicated_expected_n", "duplicated_expected_pct", "missing_n", "missing_pct",
    "consistent_n", "consistent_pct", "consistent_partial_n", "consistent_partial_pct",
    "consistent_fragmented_n", "consistent_fragmented_pct",
    "inconsistent_n", "inconsistent_pct", "inconsistent_partial_n",
    "inconsistent_partial_pct", "inconsistent_fragmented_n",
    "inconsistent_fragmented_pct", "contaminant_n", "contaminant_pct",
    "contaminant_partial_n", "contaminant_partial_pct",
    "contaminant_fragmented_n", "contaminant_fragmented_pct",
    "unknown_n", "unknown_pct", "main_species_associated_n",
    "main_species_associated_pct", "sum_file", "parser_warnings"
]

NUMERIC_COLUMNS = [c for c in COLUMNS if c.endswith("_n") or c.endswith("_pct")]
NUMERIC_COLUMNS += ["conserved_hogs", "proteins_total"]


def match_groups(pattern, text, cast=int):
    m = re.search(pattern, text, re.MULTILINE)
    return [cast(x) for x in m.groups()] if m else None


def parse_tax(path):
    result = {"selected_clade": "", "closest_species": ""}
    if not path.exists():
        return result
    section = None
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if line == ">Sampled":
            section = "selected_clade"
        elif line == ">Closest":
            section = "closest_species"
        elif line and section:
            result[section] = line
            section = None
    return result


def parse_status_table(path):
    statuses = {}
    if not path or not path.exists():
        return statuses
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("sample"):
                statuses[row["sample"]] = row.get("status", "")
    return statuses


def parse_sum(sum_path, statuses):
    sample = sum_path.stem
    text = sum_path.read_text(errors="replace")
    tax = parse_tax(sum_path.with_suffix(".tax"))
    warnings = []
    row = {c: "" for c in COLUMNS}
    row.update({
        "sample": sample,
        "status": statuses.get(sample, "COMPLETE"),
        "selected_clade": tax["selected_clade"],
        "closest_species": tax["closest_species"],
        "sum_file": str(sum_path),
    })

    m = re.search(r"^#The selected clade was\s+(.+)$", text, re.MULTILINE)
    if m and not row["selected_clade"]:
        row["selected_clade"] = m.group(1).strip()

    m = re.search(r"^#Number of conserved HOGs is:\s*(\d+)$", text, re.MULTILINE)
    if m:
        row["conserved_hogs"] = int(m.group(1))
    else:
        warnings.append("missing_conserved_hogs")

    counts = match_groups(
        r"^S:(\d+),D:(\d+)\[U:(\d+),E:(\d+)\],M:(\d+)\s*$", text
    )
    percentages = match_groups(
        r"^S:([\d.]+)%,D:([\d.]+)%\[U:([\d.]+)%,E:([\d.]+)%\],M:([\d.]+)%\s*$",
        text, float
    )
    if counts:
        s, d, du, de, missing = counts
        row.update({
            "single_n": s, "duplicated_n": d, "duplicated_unexpected_n": du,
            "duplicated_expected_n": de, "missing_n": missing,
            "complete_n": s + d,
        })
        if d != du + de:
            warnings.append("duplicated_count_mismatch")
        if row["conserved_hogs"] != "" and s + d + missing != row["conserved_hogs"]:
            warnings.append("hog_total_mismatch")
    else:
        warnings.append("missing_completeness_counts")

    if percentages:
        s, d, du, de, missing = percentages
        row.update({
            "single_pct": s, "duplicated_pct": d,
            "duplicated_unexpected_pct": du, "duplicated_expected_pct": de,
            "missing_pct": missing, "complete_pct": round(s + d, 2),
        })
    else:
        warnings.append("missing_completeness_percentages")

    m = re.search(r"^#On the whole proteome, there are\s+(\d+)\s+proteins\s*$", text, re.MULTILINE)
    if m:
        row["proteins_total"] = int(m.group(1))
    else:
        warnings.append("missing_protein_total")

    quality_counts = match_groups(
        r"^A:(\d+)\[P:(\d+),F:(\d+)\],I:(\d+)\[P:(\d+),F:(\d+)\],"
        r"C:(\d+)\[P:(\d+),F:(\d+)\],U:(\d+)\s*$", text
    )
    quality_pct = match_groups(
        r"^A:([\d.]+)%\[P:([\d.]+)%,F:([\d.]+)%\],I:([\d.]+)%"
        r"\[P:([\d.]+)%,F:([\d.]+)%\],C:([\d.]+)%"
        r"\[P:([\d.]+)%,F:([\d.]+)%\],U:([\d.]+)%\s*$", text, float
    )
    quality_names = [
        "consistent", "consistent_partial", "consistent_fragmented",
        "inconsistent", "inconsistent_partial", "inconsistent_fragmented",
        "contaminant", "contaminant_partial", "contaminant_fragmented", "unknown"
    ]
    if quality_counts:
        for name, value in zip(quality_names, quality_counts):
            row[f"{name}_n"] = value
        top_total = quality_counts[0] + quality_counts[3] + quality_counts[6] + quality_counts[9]
        if row["proteins_total"] != "" and top_total != row["proteins_total"]:
            warnings.append("protein_category_total_mismatch")
    else:
        warnings.append("missing_quality_counts")

    if quality_pct:
        for name, value in zip(quality_names, quality_pct):
            row[f"{name}_pct"] = value
    else:
        warnings.append("missing_quality_percentages")

    species_rows = []
    for line in text.splitlines():
        m = re.match(r"^(.+?)\s+(\d+)\s+(\d+)\s+([\d.]+)%\s*$", line.strip())
        if m and not line.lstrip().startswith("#"):
            species_rows.append({
                "name": m.group(1).strip(), "taxid": int(m.group(2)),
                "proteins": int(m.group(3)), "pct": float(m.group(4))
            })
    selected_species = None
    if row["closest_species"]:
        selected_species = next(
            (x for x in species_rows if x["name"] == row["closest_species"]), None
        )
    if selected_species is None and species_rows:
        selected_species = species_rows[0]
    if selected_species:
        row["closest_taxid"] = selected_species["taxid"]
        row["main_species_associated_n"] = selected_species["proteins"]
        row["main_species_associated_pct"] = selected_species["pct"]
        if not row["closest_species"]:
            row["closest_species"] = selected_species["name"]
    else:
        warnings.append("species_composition_not_parsed")

    row["parser_warnings"] = ";".join(warnings)
    return row


def write_table(path, rows, delimiter):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def numeric_summary(rows):
    output = []
    for column in NUMERIC_COLUMNS:
        values = [float(r[column]) for r in rows if r.get(column) not in ("", None)]
        if not values:
            continue
        values.sort()
        q = statistics.quantiles(values, n=4, method="inclusive") if len(values) > 1 else [values[0]] * 3
        output.append({
            "metric": column, "n": len(values), "min": min(values), "q1": q[0],
            "median": statistics.median(values), "mean": statistics.mean(values),
            "q3": q[2], "max": max(values)
        })
    return output


def main():
    parser = argparse.ArgumentParser(description="Collect OMArk .sum and .tax files into one table")
    parser.add_argument("--root", required=True, type=Path, help="Directory containing one OMArk folder per sample")
    parser.add_argument("--status", type=Path, help="Optional batch_status.tsv")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Output path prefix without extension")
    parser.add_argument("--expected", type=int, default=212, help="Expected number of samples")
    args = parser.parse_args()

    sum_files = sorted(args.root.glob("*/*.sum"))
    statuses = parse_status_table(args.status)
    rows = [parse_sum(path, statuses) for path in sum_files]
    rows.sort(key=lambda x: x["sample"].lower())

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = args.out_prefix.with_suffix(".tsv")
    csv_path = args.out_prefix.with_suffix(".csv")
    stats_path = args.out_prefix.parent / f"{args.out_prefix.name}_descriptive_stats.tsv"
    report_path = args.out_prefix.parent / f"{args.out_prefix.name}_validation.txt"

    write_table(tsv_path, rows, "\t")
    write_table(csv_path, rows, ",")

    stats = numeric_summary(rows)
    with stats_path.open("w", newline="") as handle:
        fields = ["metric", "n", "min", "q1", "median", "mean", "q3", "max"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(stats)

    warning_rows = [r for r in rows if r["parser_warnings"]]
    incomplete_rows = [r for r in rows if r["status"] != "COMPLETE"]
    with report_path.open("w") as handle:
        handle.write(f"Expected samples: {args.expected}\n")
        handle.write(f"SUM files found: {len(sum_files)}\n")
        handle.write(f"Rows written: {len(rows)}\n")
        handle.write(f"Rows with parser warnings: {len(warning_rows)}\n")
        handle.write(f"Rows not marked COMPLETE: {len(incomplete_rows)}\n")
        handle.write(f"Count check: {'PASS' if len(rows) == args.expected else 'FAIL'}\n")
        if warning_rows:
            handle.write("\nParser warnings:\n")
            for row in warning_rows:
                handle.write(f"{row['sample']}\t{row['parser_warnings']}\n")
        if incomplete_rows:
            handle.write("\nNon-complete statuses:\n")
            for row in incomplete_rows:
                handle.write(f"{row['sample']}\t{row['status']}\n")

    print(f"Rows: {len(rows)}")
    print(f"Main TSV: {tsv_path}")
    print(f"CSV: {csv_path}")
    print(f"Descriptive statistics: {stats_path}")
    print(f"Validation report: {report_path}")
    if len(rows) != args.expected or warning_rows or incomplete_rows:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
