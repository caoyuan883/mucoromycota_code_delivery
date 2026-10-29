#!/usr/bin/env python3

import csv
import argparse
from pathlib import Path


def read_fasta_lengths_and_gc(fp):
    lengths = []
    total_len = 0
    gc = 0
    seq_len = 0

    with open(fp, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line:
                continue
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


def parse_sample_name(fp):
    sample = fp.name
    if sample.endswith(".fna"):
        sample = sample[:-4]
    if sample.endswith("_genomic"):
        sample = sample[:-8]
    return sample


def collect_assembly_stats(genome_dir):
    rows = []
    fasta_files = sorted(genome_dir.glob("*.fna"))

    for fp in fasta_files:
        sample = parse_sample_name(fp)
        lengths, total_len, gc_percent = read_fasta_lengths_and_gc(fp)
        n50, l50 = nx_lx(lengths, 0.5)
        n90, l90 = nx_lx(lengths, 0.9)

        rows.append({
            "sample": sample,
            "genome_size_bp": total_len,
            "genome_size_mb": round(total_len / 1_000_000, 4),
            "GC_percent": gc_percent,
            "scaffold_count": len(lengths),
            "longest_scaffold_bp": max(lengths) if lengths else 0,
            "N50_bp": n50,
            "L50": l50,
            "N90_bp": n90,
            "L90": l90,
        })

    return rows


def read_csv_rows(csv_file):
    rows = []
    with open(csv_file, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows


def write_csv_rows(csv_file, rows, fieldnames):
    with open(csv_file, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Build assembly statistics and merge them with a BUSCO summary table."
    )
    parser.add_argument(
        "--genome-dir",
        required=True,
        help="Directory containing genome FASTA files (*.fna)"
    )
    parser.add_argument(
        "--busco-csv",
        required=True,
        help="BUSCO summary table in wide CSV format"
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for generated CSV files"
    )
    args = parser.parse_args()

    genome_dir = Path(args.genome_dir).expanduser().resolve()
    busco_csv = Path(args.busco_csv).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()

    if not genome_dir.is_dir():
        raise SystemExit(f"ERROR: genome directory not found: {genome_dir}")
    if not busco_csv.is_file():
        raise SystemExit(f"ERROR: BUSCO summary CSV not found: {busco_csv}")

    out_dir.mkdir(parents=True, exist_ok=True)

    assembly_csv = out_dir / "genome_basic_stats.csv"
    master_csv = out_dir / "genome_qc_master_table.csv"

    assembly_rows = collect_assembly_stats(genome_dir)

    assembly_fieldnames = [
        "sample", "genome_size_bp", "genome_size_mb", "GC_percent",
        "scaffold_count", "longest_scaffold_bp", "N50_bp", "L50", "N90_bp", "L90"
    ]
    write_csv_rows(assembly_csv, assembly_rows, assembly_fieldnames)

    busco_rows = read_csv_rows(busco_csv)
    asm_map = {r["sample"]: r for r in assembly_rows}

    merged = []
    for r in busco_rows:
        sample = r["sample"]
        m = dict(r)
        if sample in asm_map:
            m.update(asm_map[sample])
        merged.append(m)

    busco_fields = list(busco_rows[0].keys()) if busco_rows else ["sample"]
    asm_fields = [
        "genome_size_bp", "genome_size_mb", "GC_percent",
        "scaffold_count", "longest_scaffold_bp", "N50_bp", "L50", "N90_bp", "L90"
    ]
    fieldnames = busco_fields + [x for x in asm_fields if x not in busco_fields]

    write_csv_rows(master_csv, merged, fieldnames)

    print(f"written: {assembly_csv}")
    print(f"written: {master_csv}")


if __name__ == "__main__":
    main()