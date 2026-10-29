#!/usr/bin/env python3
"""
Filters run-level SRA metadata to retain only paired-end Illumina RNA-seq
libraries, keeping up to six runs per species (in the order returned by
Entrez), following the selection criteria described in Nenasheva et al.
(2025, Scientific Data) for diatom genome annotation.
"""

import pandas as pd

INPUT_FILE = "sra_detailed_runinfo.csv"
OUTPUT_FILE = "sci_grade_download_list_v2.csv"
MAX_RUNS_PER_SPECIES = 6

df = pd.read_csv(INPUT_FILE)
df["AvgLength"] = pd.to_numeric(df["AvgLength"], errors="coerce")

filtered = df[df["LibraryLayout"] == "PAIRED"].copy()
final_runs = filtered.groupby("Species").head(MAX_RUNS_PER_SPECIES)

final_runs.to_csv(OUTPUT_FILE, index=False)

n_species = final_runs["Species"].nunique()
n_runs = len(final_runs)

print("=" * 50)
print("RNA-seq run filtering report")
print("=" * 50)
print(f"Species with RNA-seq (BRAKER3 candidates): {n_species}")
print(f"Total runs in final download list: {n_runs}")
if n_species > 0:
    print(f"Average runs per species: {n_runs / n_species:.2f}")
else:
    print("Average runs per species: 0")
print("=" * 50)
print("Note: co-culture libraries are not excluded (following reference method).")
print("Note: alignment-rate filtering is applied downstream, after HISAT2 mapping.")
print("Note: no minimum read-length threshold is applied.")
