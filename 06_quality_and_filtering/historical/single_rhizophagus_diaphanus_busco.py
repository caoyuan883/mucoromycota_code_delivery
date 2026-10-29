import csv, re
from pathlib import Path

FNA = Path("/home/py/Desktop/12_JGI/Rhizophagus_diaphanus_genomic.fna")
PROJECT = Path("/home/py/bioinfo/mucoromycota_busco_local_13_JGI_20260627_143655")
OUTDIR = PROJECT / "rhizophagus_diaphanus_single_extract"
OUTDIR.mkdir(parents=True, exist_ok=True)

TARGET = "Rhizophagus_diaphanus"

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

def parse_busco_txt(txt_path):
    text = Path(txt_path).read_text(encoding="utf-8", errors="ignore")
    out = {
        "dataset": "",
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

    pats = {
        "C_n": r"Complete BUSCOs \(C\)\s+(\d+)",
        "S_n": r"Complete and single-copy BUSCOs \(S\)\s+(\d+)",
        "D_n": r"Complete and duplicated BUSCOs \(D\)\s+(\d+)",
        "F_n": r"Fragmented BUSCOs \(F\)\s+(\d+)",
        "M_n": r"Missing BUSCOs \(M\)\s+(\d+)",
        "n_buscos_alt": r"Total BUSCO groups searched\s+(\d+)",
        "busco_version": r"BUSCO version is:\s*([^\n]+)",
        "dataset": r"The lineage dataset is:\s*([^\n]+)",
    }
    for k, p in pats.items():
        mm = re.search(p, text)
        if mm:
            if k == "n_buscos_alt" and not out["n_buscos"]:
                out["n_buscos"] = mm.group(1)
            elif k in out:
                out[k] = mm.group(1).strip()
    return out

def pick_summary_txt(sample_dir):
    txts = sorted(sample_dir.rglob("short_summary*.txt"))
    for x in txts:
        if "specific" in x.name:
            return x
    return txts[0] if txts else None

rows = []
for lineage_group in ["results_mucoromycota", "results_fungi"]:
    root = PROJECT / lineage_group
    if not root.exists():
        continue
    candidates = [d for d in root.iterdir() if d.is_dir() and TARGET.lower() in d.name.lower()]
    for sample_dir in candidates:
        txt = pick_summary_txt(sample_dir)
        fulls = sorted(sample_dir.rglob("full_table.tsv"))
        miss = sorted(sample_dir.rglob("missing_busco_list.tsv"))
        row = {
            "target_name": TARGET,
            "matched_sample_dir": sample_dir.name,
            "lineage_group": lineage_group,
            "summary_txt": str(txt) if txt else "",
            "full_table_tsv": str(fulls[0]) if fulls else "",
            "missing_busco_list": str(miss[0]) if miss else "",
            "status": "ok" if txt else "missing_summary_txt"
        }
        if txt:
            row.update(parse_busco_txt(txt))
        rows.append(row)

lengths, total_len, gc_percent = read_fasta_lengths_and_gc(FNA)
n50, l50 = nx_lx(lengths, 0.5)
n90, l90 = nx_lx(lengths, 0.9)

asm_row = {
    "sample": TARGET,
    "fna_path": str(FNA),
    "genome_size_bp": total_len,
    "genome_size_mb": round(total_len / 1_000_000, 4),
    "GC_percent": gc_percent,
    "scaffold_count": len(lengths),
    "longest_scaffold_bp": max(lengths) if lengths else 0,
    "N50_bp": n50,
    "L50": l50,
    "N90_bp": n90,
    "L90": l90,
}

with open(OUTDIR / "01_single_busco_Rhizophagus_diaphanus.csv", "w", newline="", encoding="utf-8") as out:
    fieldnames = [
        "target_name","matched_sample_dir","lineage_group","dataset","busco_version",
        "C_percent","S_percent","D_percent","F_percent","M_percent",
        "C_n","S_n","D_n","F_n","M_n","n_buscos",
        "summary_txt","full_table_tsv","missing_busco_list","status"
    ]
    w = csv.DictWriter(out, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

with open(OUTDIR / "02_single_genome_qc_Rhizophagus_diaphanus.csv", "w", newline="", encoding="utf-8") as out:
    fieldnames = [
        "sample","fna_path","genome_size_bp","genome_size_mb","GC_percent",
        "scaffold_count","longest_scaffold_bp","N50_bp","L50","N90_bp","L90"
    ]
    w = csv.DictWriter(out, fieldnames=fieldnames)
    w.writeheader()
    w.writerow(asm_row)

print("written:", OUTDIR / "01_single_busco_Rhizophagus_diaphanus.csv")
print("written:", OUTDIR / "02_single_genome_qc_Rhizophagus_diaphanus.csv")
print("matched_rows:", len(rows))
