#!/usr/bin/env python3
"""
Aggregates the run-level RNA-seq download table (with genome paths) into a
species-level processing plan: one row per species, listing all associated
SRA run accessions, library layouts, and the matched genome path.
"""

import argparse
from collections import defaultdict

import pandas as pd

DEFAULT_INPUT = "18_tables_sra_run_download_list_with_genome_path.xlsx"
DEFAULT_OUTPUT = "19_tables_species_level_rnaseq_plan.xlsx"


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate run-level table into species-level RNA-seq plan."
    )
    parser.add_argument("--input", default=DEFAULT_INPUT,
                         help="Input xlsx table (run-level, with GenomePath)")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                         help="Output xlsx table (species-level plan)")
    args = parser.parse_args()

    df = pd.read_excel(args.input)

    species_info = defaultdict(lambda: {"runs": [], "layouts": set(), "genome": None})

    for _, row in df.iterrows():
        sp = str(row["Species"]).strip()
        species_info[sp]["runs"].append(str(row["Run_Accession"]).strip())
        species_info[sp]["layouts"].add(str(row["LibraryLayout"]).strip())
        species_info[sp]["genome"] = str(row["GenomePath"]).strip()

    records = []
    for sp in sorted(species_info):
        info = species_info[sp]
        records.append({
            "Species": sp,
            "RunCount": len(info["runs"]),
            "Runs": ";".join(info["runs"]),
            "LibraryLayout": ";".join(sorted(info["layouts"])),
            "GenomePath": info["genome"],
        })

    out_df = pd.DataFrame(records)
    out_df.to_excel(args.output, index=False)

    print(f"Output table: {args.output}")
    print(f"Species count: {len(species_info)}")
    print(f"Total run records: {sum(len(v['runs']) for v in species_info.values())}")


if __name__ == "__main__":
    main()
