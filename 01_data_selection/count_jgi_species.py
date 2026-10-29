#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script : count_jgi_species.py
Purpose: Normalize species names from JGI metadata and generate a deduplicated 
         count of unique species/strains.
"""

import csv
import re
import sys
import argparse
from collections import Counter
from pathlib import Path

def normalize_name(name):
    """
    Normalizes the taxonomic name to 'Genus species' or 'Genus sp. identifier'.
    """
    s = str(name).strip()
    s = re.sub(r'\s+', ' ', s)
    if not s:
        return None

    parts = s.split(' ')
    if len(parts) < 2:
        return s

    genus = parts[0]
    second = parts[1]

    # Preserve strain-specific identifiers if "sp." is present
    if second == "sp.":
        if len(parts) >= 3:
            return f"{genus} sp. {parts[2]}"
        else:
            return f"{genus} sp."
    else:
        return f"{genus} {second}"

def main():
    parser = argparse.ArgumentParser(description="Normalize and count unique species from JGI dataset.")
    parser.add_argument("--input-tsv", required=True, help="Path to the raw JGI TSV/CSV file.")
    parser.add_argument("--out-tsv", required=True, help="Path to the output species count TSV file.")
    args = parser.parse_args()

    input_file = Path(args.input_tsv).resolve()
    output_file = Path(args.out_tsv).resolve()

    if not input_file.is_file():
        sys.exit(f"ERROR: Input file not found: {input_file}")

    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)

    counts = Counter()
    raw_total = 0

    # Read input file, inferring delimiter (tab or comma)
    with open(input_file, "r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
            delim = dialect.delimiter
        except:
            delim = "\t"

        reader = csv.reader(f, delimiter=delim)
        for row in reader:
            if not row:
                continue
            first_col = row[0].strip()
            if not first_col:
                continue
            
            raw_total += 1
            species = normalize_name(first_col)
            if species:
                counts[species] += 1

    # Write output
    with open(output_file, "w", encoding="utf-8", newline="") as out:
        out.write("species_name\tcount\n")
        # Sort primarily by descending count, secondarily by ascending name
        for species_name, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
            out.write(f"{species_name}\t{count}\n")

    print("===== JGI Species Count Summary =====")
    print(f"Total rows counted : {raw_total}")
    print(f"Unique species     : {len(counts)}")
    print(f"Output written to  : {output_file}")

if __name__ == "__main__":
    main()
