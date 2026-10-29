import json
import csv
import re
from pathlib import Path

base = Path.home() / "bioinfo" / "protein_busco_161"
dirs = {
    "mucoromycota_odb12": base / "results_mucoromycota",
    "fungi_odb12": base / "results_fungi",
}

out_long = base / "protein_busco_summary_long.csv"
out_wide = base / "protein_busco_summary_wide.csv"
out_issues = base / "protein_busco_parse_issues.tsv"
out_paths = base / "protein_busco_result_paths.tsv"

all_samples = set()
results_dict = {}
paths_rows = []
issue_rows = []

def norm_value(v):
    if v is None:
        return ""
    return str(v)

def parse_json_summary(json_file):
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)

    results = data.get("results", {}).get("Results", {})
    if not results:
        return None

    return {
        "C_percent": norm_value(results.get("Complete BUSCOs %", "")),
        "S_percent": norm_value(results.get("Single-copy BUSCOs %", "")),
        "D_percent": norm_value(results.get("Multi-copy BUSCOs %", "")),
        "F_percent": norm_value(results.get("Fragmented BUSCOs %", "")),
        "M_percent": norm_value(results.get("Missing BUSCOs %", "")),
        "n_buscos": norm_value(results.get("Total BUSCO groups searched", "")),
        "source_type": "json",
        "source_file": str(json_file),
    }

def parse_txt_summary(txt_file):
    text = txt_file.read_text(encoding="utf-8", errors="ignore")
    target_line = None
    for line in text.splitlines():
        if all(x in line for x in ["C:", "S:", "D:", "F:", "M:", "n:"]):
            target_line = line.strip()
            break
    if target_line is None:
        return None

    def grab(pattern):
        m = re.search(pattern, target_line)
        return m.group(1) if m else ""

    return {
        "C_percent": grab(r"C:([0-9.]+)%"),
        "S_percent": grab(r"S:([0-9.]+)%"),
        "D_percent": grab(r"D:([0-9.]+)%"),
        "F_percent": grab(r"F:([0-9.]+)%"),
        "M_percent": grab(r"M:([0-9.]+)%"),
        "n_buscos": grab(r"n:([0-9]+)"),
        "source_type": "txt",
        "source_file": str(txt_file),
    }

for lineage, root in dirs.items():
    if not root.exists():
        issue_rows.append({
            "sample": "",
            "lineage": lineage,
            "issue_type": "ROOT_DIR_MISSING",
            "details": str(root)
        })
        continue

    for sample_dir in sorted([x for x in root.iterdir() if x.is_dir()]):
        sample = sample_dir.name
        all_samples.add(sample)
        results_dict.setdefault(sample, {})

        json_files = sorted(sample_dir.rglob("short_summary*.json"))
        txt_files = sorted(sample_dir.rglob("short_summary*.txt"))
        full_tables = sorted(sample_dir.rglob("full_table.tsv"))
        missing_lists = sorted(sample_dir.rglob("missing_busco_list.tsv"))

        paths_rows.append({
            "sample": sample,
            "lineage": lineage,
            "sample_dir": str(sample_dir),
            "json_summary_file": str(json_files[0]) if json_files else "",
            "txt_summary_file": str(txt_files[0]) if txt_files else "",
            "full_table_tsv": str(full_tables[0]) if full_tables else "",
            "missing_busco_list_tsv": str(missing_lists[0]) if missing_lists else "",
        })

        parsed = None

        if json_files:
            try:
                parsed = parse_json_summary(json_files[0])
            except Exception as e:
                issue_rows.append({
                    "sample": sample,
                    "lineage": lineage,
                    "issue_type": "JSON_PARSE_FAILED",
                    "details": f"{json_files[0]} | {e}"
                })

        if parsed is None and txt_files:
            try:
                parsed = parse_txt_summary(txt_files[0])
            except Exception as e:
                issue_rows.append({
                    "sample": sample,
                    "lineage": lineage,
                    "issue_type": "TXT_PARSE_FAILED",
                    "details": f"{txt_files[0]} | {e}"
                })

        if parsed is None:
            issue_rows.append({
                "sample": sample,
                "lineage": lineage,
                "issue_type": "NO_PARSEABLE_SUMMARY",
                "details": str(sample_dir)
            })
            continue

        results_dict[sample][lineage] = parsed

with open(out_long, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow([
        "sample", "lineage",
        "C_percent", "S_percent", "D_percent", "F_percent", "M_percent", "n_buscos",
        "source_type", "source_file"
    ])
    for sample in sorted(all_samples):
        for lineage in ["mucoromycota_odb12", "fungi_odb12"]:
            r = results_dict.get(sample, {}).get(lineage, {})
            writer.writerow([
                sample, lineage,
                r.get("C_percent", ""),
                r.get("S_percent", ""),
                r.get("D_percent", ""),
                r.get("F_percent", ""),
                r.get("M_percent", ""),
                r.get("n_buscos", ""),
                r.get("source_type", ""),
                r.get("source_file", "")
            ])

with open(out_wide, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow([
        "sample",
        "C_mucoromycota", "S_mucoromycota", "D_mucoromycota", "F_mucoromycota", "M_mucoromycota", "n_mucoromycota",
        "C_fungi", "S_fungi", "D_fungi", "F_fungi", "M_fungi", "n_fungi"
    ])
    for sample in sorted(all_samples):
        a = results_dict.get(sample, {}).get("mucoromycota_odb12", {})
        b = results_dict.get(sample, {}).get("fungi_odb12", {})
        writer.writerow([
            sample,
            a.get("C_percent", ""), a.get("S_percent", ""), a.get("D_percent", ""), a.get("F_percent", ""), a.get("M_percent", ""), a.get("n_buscos", ""),
            b.get("C_percent", ""), b.get("S_percent", ""), b.get("D_percent", ""), b.get("F_percent", ""), b.get("M_percent", ""), b.get("n_buscos", "")
        ])

with open(out_issues, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=["sample", "lineage", "issue_type", "details"], delimiter="\t")
    writer.writeheader()
    for row in issue_rows:
        writer.writerow(row)

with open(out_paths, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "sample", "lineage", "sample_dir",
            "json_summary_file", "txt_summary_file",
            "full_table_tsv", "missing_busco_list_tsv"
        ],
        delimiter="\t"
    )
    writer.writeheader()
    for row in sorted(paths_rows, key=lambda x: (x["sample"], x["lineage"])):
        writer.writerow(row)

print(f"written: {out_long}")
print(f"written: {out_wide}")
print(f"written: {out_issues}")
print(f"written: {out_paths}")
print(f"total_samples: {len(all_samples)}")
print(f"parsed_mucoromycota: {sum(1 for s in all_samples if 'mucoromycota_odb12' in results_dict.get(s, {}))}")
print(f"parsed_fungi: {sum(1 for s in all_samples if 'fungi_odb12' in results_dict.get(s, {}))}")
print(f"issues: {len(issue_rows)}")
