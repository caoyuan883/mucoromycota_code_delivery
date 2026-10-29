#!/usr/bin/env python3
"""
Maps each species in the RNA-seq run download list to its corresponding
local genome FASTA file path, based on filename matching
(Genus_species_genomic.fna). Produces an augmented table with an added
GenomePath column, plus a list of species with no matching genome file.
"""

import argparse
from pathlib import Path

import pandas as pd

DEFAULT_INPUT = "17_tables_sra_run_download_list.xlsx"
DEFAULT_GENOME_DIR = "genomes"
DEFAULT_OUTPUT = "18_tables_sra_run_download_list_with_genome_path.xlsx"
DEFAULT_UNMATCHED = "18_unmatched_species.txt"


def build_genome_map(genome_dir: Path) -> dict:
    genome_map = {}
    for f in genome_dir.glob("*"):
        if f.is_file():
            genome_map[f.name] = str(f.resolve())
    return genome_map


def main():
    parser = argparse.ArgumentParser(
        description="Map species to local genome file paths."
    )
    parser.add_argument("--input", default=DEFAULT_INPUT,
                         help="Input xlsx table with a 'Species' column")
    parser.add_argument("--genome-dir", default=DEFAULT_GENOME_DIR,
                         help="Directory containing genome FASTA files")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                         help="Output xlsx table with GenomePath column added")
    parser.add_argument("--unmatched", default=DEFAULT_UNMATCHED,
                         help="Output text file listing unmatched species")
    args = parser.parse_args()

    csv_in = Path(args.input)
    genome_dir = Path(args.genome_dir)
    xlsx_out = Path(args.output)
    unmatched_out = Path(args.unmatched)

    df = pd.read_excel(csv_in)
    genome_map = build_genome_map(genome_dir)

    matched = 0
    unmatched = []
    genome_paths = []

    for species in df["Species"].astype(str):
        species_key = species.strip().replace(" ", "_")
        target_name = f"{species_key}_genomic.fna"
        genome_path = genome_map.get(target_name, "")
        if genome_path:
            matched += 1
        else:
            unmatched.append(species)
        genome_paths.append(genome_path)

    df["GenomePath"] = genome_paths
    df.to_excel(xlsx_out, index=False)

    with open(unmatched_out, "w", encoding="utf-8") as f:
        for s in unmatched:
            f.write(s + "\n")

    print(f"Output table: {xlsx_out}")
    print(f"Matched: {matched}")
    print(f"Unmatched: {len(unmatched)}")
    print(f"Unmatched list: {unmatched_out}")


if __name__ == "__main__":
    main()
