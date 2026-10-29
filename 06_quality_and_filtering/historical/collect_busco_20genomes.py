import os, csv, json, re
from pathlib import Path

BASE = Path("/home/py/bioinfo")
PROJECT_DIRS = [
    "mucoromycota_busco_local_7genomic_20260625_033828",
    "mucoromycota_busco_local_13_JGI_20260627_143655",
    "mucoromycota_busco_local_Benjaminiella_poitrasii_JGI_20260630_021824",
]

FNA_SOURCES = [
    Path("/home/py/Desktop/7_genomic"),
    Path("/home/py/Desktop/12_JGI"),
    Path("/home/py/Desktop/Benjaminiella_poitrasii_JGI_genomic.fna"),
]

OUT_DIR = BASE / "busco_manual_merge_20genomes_fix2_txt_first"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def read_fasta_lengths_and_gc(fp):
    lengths = []
    total_len = 0
    gc = 0
    seq_len = 0
    with open(fp, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                    seq_len = 0
            else:
                s = line.strip().upper()
                if not s:
                    continue
                seq_len += len(s)
                total_len += len(s)
                gc += s.count("G") + s.count("C")
        if seq_len > 0:
            lengths.append(seq_len)
    gc_percent = round((gc / total_len * 100), 4) if total_len > 0 else 0
    return lengths, total_len, gc_percent

def nx_lx(lengths, frac):
    if not lengths:
        return 0, 0
    arr = sorted(lengths, reverse=True)
    target = sum(arr) * frac
    csum = 0
    for i, x in enumerate(arr, start=1):
        csum += x
        if csum >= target:
            return x, i
    return 0, 0

def build_fna_index():
    idx = {}
    for src in FNA_SOURCES:
        if src.is_dir():
            for fna in src.glob("*.fna"):
                sample = fna.name
                if sample.endswith(".fna"):
                    sample = sample[:-4]
                if sample.endswith("_genomic"):
                    sample = sample[:-8]
                idx[sample] = str(fna)
        elif src.is_file() and src.name.endswith(".fna"):
            sample = src.name
            if sample.endswith(".fna"):
                sample = sample[:-4]
            if sample.endswith("_genomic"):
                sample = sample[:-8]
            idx[sample] = str(src)
    return idx

def parse_busco_txt(txt_path):
    text = Path(txt_path).read_text(encoding="utf-8", errors="ignore")

    out = {
        "dataset": "",
        "parent_dataset": "",
        "busco_version": "",
        "C_percent": "", "S_percent": "", "D_percent": "", "F_percent": "", "M_percent": "",
        "C_n": "", "S_n": "", "D_n": "", "F_n": "", "M_n": "", "n_buscos": "",
    }

    m = re.search(r"C:(\d+(?:\.\d+)?)%\[S:(\d+(?:\.\d+)?)%,D:(\d+(?:\.\d+)?)%\],F:(\d+(?:\.\d+)?)%,M:(\d+(?:\.\d+)?)%,n:(\d+)", text)
    if m:
        out["C_percent"] = m.group(1)
        out["S_percent"] = m.group(2)
        out["D_percent"] = m.group(3)
        out["F_percent"] = m.group(4)
        out["M_percent"] = m.group(5)
        out["n_buscos"] = m.group(6)

    patterns = {
        "C_n": r"Complete BUSCOs \(C\)\s+(\d+)",
        "S_n": r"Complete and single-copy BUSCOs \(S\)\s+(\d+)",
        "D_n": r"Complete and duplicated BUSCOs \(D\)\s+(\d+)",
        "F_n": r"Fragmented BUSCOs \(F\)\s+(\d+)",
        "M_n": r"Missing BUSCOs \(M\)\s+(\d+)",
        "n_buscos_alt": r"Total BUSCO groups searched\s+(\d+)",
        "busco_version": r"BUSCO version is:\s*([^\n]+)",
        "dataset": r"The lineage dataset is:\s*([^\n]+)",
        "parent_dataset": r"Placement file versions:\s*([^\n]+)",
    }

    for key, pat in patterns.items():
        mm = re.search(pat, text)
        if mm:
            if key == "n_buscos_alt" and not out["n_buscos"]:
                out["n_buscos"] = mm.group(1)
            elif key in out:
                out[key] = mm.group(1).strip()

    filled = any(out[k] != "" for k in ["C_percent","S_percent","D_percent","F_percent","M_percent","n_buscos"])
    return out, filled, text[:500].replace("\n", " | ")

def get_first(d, keys, default=""):
    for k in keys:
        if isinstance(d, dict) and k in d and d[k] not in [None, ""]:
            return d[k]
    return default

def parse_busco_json(js_path):
    with open(js_path, "r", encoding="utf-8", errors="ignore") as f:
        data = json.load(f)

    res = {}
    if isinstance(data.get("results"), dict):
        res = data["results"].get("Results", {}) or data["results"].get("results", {}) or {}
    versions = data.get("versions", {}) if isinstance(data.get("versions"), dict) else {}

    out = {
        "dataset": get_first(data, ["lineage_dataset", "dataset", "odb_dataset"], ""),
        "parent_dataset": str(get_first(data, ["parent_dataset"], "")),
        "busco_version": get_first(versions, ["busco"], get_first(data, ["busco_version", "version"], "")),
        "C_percent": get_first(res, ["Complete BUSCOs %", "C"]),
        "S_percent": get_first(res, ["Single-copy BUSCOs %", "S"]),
        "D_percent": get_first(res, ["Multi-copy BUSCOs %", "Duplicated BUSCOs %", "D"]),
        "F_percent": get_first(res, ["Fragmented BUSCOs %", "F"]),
        "M_percent": get_first(res, ["Missing BUSCOs %", "M"]),
        "C_n": get_first(res, ["Complete BUSCOs", "Complete"]),
        "S_n": get_first(res, ["Single-copy BUSCOs", "Single copy BUSCOs"]),
        "D_n": get_first(res, ["Multi-copy BUSCOs", "Duplicated BUSCOs"]),
        "F_n": get_first(res, ["Fragmented BUSCOs", "Fragmented"]),
        "M_n": get_first(res, ["Missing BUSCOs", "Missing"]),
        "n_buscos": get_first(res, ["Total BUSCO groups searched", "n", "Total"]),
    }
    filled = any(out[k] != "" for k in ["C_percent","S_percent","D_percent","F_percent","M_percent","n_buscos"])
    return out, filled

fna_index = build_fna_index()

busco_rows = []
path_rows = []
asm_rows = []
asm_done = set()
debug_rows = []

for proj in PROJECT_DIRS:
    proj_path = BASE / proj
    for lineage_group in ["results_mucoromycota", "results_fungi"]:
        root = proj_path / lineage_group
        if not root.exists():
            continue

        for sample_dir in sorted([x for x in root.iterdir() if x.is_dir()]):
            sample = sample_dir.name

            json_files = sorted(sample_dir.rglob("short_summary*.json"))
            txt_files = sorted(sample_dir.rglob("short_summary*.txt"))
            full_tables = sorted(sample_dir.rglob("full_table.tsv"))
            missing_lists = sorted(sample_dir.rglob("missing_busco_list.tsv"))

            chosen_txt = None
            for tf in txt_files:
                if "specific" in tf.name:
                    chosen_txt = tf
                    break
            if chosen_txt is None and txt_files:
                chosen_txt = txt_files[0]

            chosen_json = None
            for jf in json_files:
                if "specific" in jf.name:
                    chosen_json = jf
                    break
            if chosen_json is None and json_files:
                chosen_json = json_files[0]

            row = {
                "project_folder": proj,
                "lineage_group": lineage_group,
                "sample": sample,
                "summary_json": str(chosen_json) if chosen_json else "",
                "summary_txt": str(chosen_txt) if chosen_txt else "",
                "full_table_tsv_files": ";".join(str(x) for x in full_tables),
                "missing_busco_list_files": ";".join(str(x) for x in missing_lists),
            }

            parse_source = "none"
            preview = ""

            if chosen_txt and Path(chosen_txt).exists():
                parsed_txt, ok_txt, preview = parse_busco_txt(chosen_txt)
                if ok_txt:
                    row.update(parsed_txt)
                    row["status"] = "ok"
                    parse_source = "txt"
                else:
                    row.update(parsed_txt)

            if parse_source == "none" and chosen_json and Path(chosen_json).exists():
                parsed_json, ok_json = parse_busco_json(chosen_json)
                if ok_json:
                    row.update(parsed_json)
                    row["status"] = "ok"
                    parse_source = "json"
                else:
                    row.update(parsed_json)

            if parse_source == "none":
                row.update({
                    "dataset": row.get("dataset",""),
                    "parent_dataset": row.get("parent_dataset",""),
                    "busco_version": row.get("busco_version",""),
                    "C_percent": row.get("C_percent",""),
                    "S_percent": row.get("S_percent",""),
                    "D_percent": row.get("D_percent",""),
                    "F_percent": row.get("F_percent",""),
                    "M_percent": row.get("M_percent",""),
                    "C_n": row.get("C_n",""),
                    "S_n": row.get("S_n",""),
                    "D_n": row.get("D_n",""),
                    "F_n": row.get("F_n",""),
                    "M_n": row.get("M_n",""),
                    "n_buscos": row.get("n_buscos",""),
                    "status": "parse_failed"
                })

            row["parse_source"] = parse_source

            debug_rows.append({
                "project_folder": proj,
                "lineage_group": lineage_group,
                "sample": sample,
                "summary_txt": row["summary_txt"],
                "summary_json": row["summary_json"],
                "parse_source": parse_source,
                "txt_preview": preview,
                "C_percent": row.get("C_percent",""),
                "S_percent": row.get("S_percent",""),
                "D_percent": row.get("D_percent",""),
                "F_percent": row.get("F_percent",""),
                "M_percent": row.get("M_percent",""),
                "n_buscos": row.get("n_buscos",""),
                "status": row["status"]
            })

            busco_rows.append({k: row.get(k, "") for k in [
                "project_folder","lineage_group","sample","dataset","parent_dataset","busco_version",
                "C_percent","S_percent","D_percent","F_percent","M_percent",
                "C_n","S_n","D_n","F_n","M_n","n_buscos","status","parse_source"
            ]})

            path_rows.append({k: row.get(k, "") for k in [
                "project_folder","lineage_group","sample","summary_txt","summary_json",
                "full_table_tsv_files","missing_busco_list_files","parse_source"
            ]})

            if sample not in asm_done:
                asm_done.add(sample)
                fna_path = fna_index.get(sample, "")
                if fna_path and Path(fna_path).exists():
                    lengths, total_len, gc_percent = read_fasta_lengths_and_gc(fna_path)
                    n50, l50 = nx_lx(lengths, 0.5)
                    n90, l90 = nx_lx(lengths, 0.9)
                    asm_rows.append({
                        "sample": sample,
                        "fna_path": fna_path,
                        "genome_size_bp": total_len,
                        "genome_size_mb": round(total_len / 1_000_000, 4),
                        "GC_percent": gc_percent,
                        "scaffold_count": len(lengths),
                        "longest_scaffold_bp": max(lengths) if lengths else 0,
                        "N50_bp": n50,
                        "L50": l50,
                        "N90_bp": n90,
                        "L90": l90,
                        "assembly_status": "ok"
                    })
                else:
                    asm_rows.append({
                        "sample": sample,
                        "fna_path": "",
                        "genome_size_bp": "",
                        "genome_size_mb": "",
                        "GC_percent": "",
                        "scaffold_count": "",
                        "longest_scaffold_bp": "",
                        "N50_bp": "",
                        "L50": "",
                        "N90_bp": "",
                        "L90": "",
                        "assembly_status": "fna_not_found"
                    })

asm_map = {r["sample"]: r for r in asm_rows}
master_rows = []
for r in busco_rows:
    m = dict(r)
    if r["sample"] in asm_map:
        for k, v in asm_map[r["sample"]].items():
            if k != "sample":
                m[k] = v
    master_rows.append(m)

def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

write_csv(
    OUT_DIR / "00_busco_parse_debug.csv",
    debug_rows,
    ["project_folder","lineage_group","sample","summary_txt","summary_json","parse_source",
     "txt_preview","C_percent","S_percent","D_percent","F_percent","M_percent","n_buscos","status"]
)

write_csv(
    OUT_DIR / "01_busco_summary_20genomes_wide.csv",
    busco_rows,
    ["project_folder","lineage_group","sample","dataset","parent_dataset","busco_version",
     "C_percent","S_percent","D_percent","F_percent","M_percent",
     "C_n","S_n","D_n","F_n","M_n","n_buscos","status","parse_source"]
)

write_csv(
    OUT_DIR / "02_busco_result_paths_20genomes.csv",
    path_rows,
    ["project_folder","lineage_group","sample","summary_txt","summary_json",
     "full_table_tsv_files","missing_busco_list_files","parse_source"]
)

write_csv(
    OUT_DIR / "03_genome_basic_stats_20genomes.csv",
    asm_rows,
    ["sample","fna_path","genome_size_bp","genome_size_mb","GC_percent",
     "scaffold_count","longest_scaffold_bp","N50_bp","L50","N90_bp","L90","assembly_status"]
)

write_csv(
    OUT_DIR / "04_genome_qc_master_table_20genomes.csv",
    master_rows,
    ["project_folder","lineage_group","sample","dataset","parent_dataset","busco_version",
     "C_percent","S_percent","D_percent","F_percent","M_percent",
     "C_n","S_n","D_n","F_n","M_n","n_buscos","status","parse_source",
     "fna_path","genome_size_bp","genome_size_mb","GC_percent",
     "scaffold_count","longest_scaffold_bp","N50_bp","L50","N90_bp","L90","assembly_status"]
)

print("written:", OUT_DIR / "00_busco_parse_debug.csv")
print("written:", OUT_DIR / "04_genome_qc_master_table_20genomes.csv")
print("busco_rows:", len(busco_rows))
print("unique_samples:", len(asm_rows))
