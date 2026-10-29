#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path
import re

input_file = Path.home() / "Desktop" / "01_tables_raw_NCBI_Mucoromycota.tsv"
output_dir = Path.home() / "Desktop" / "ncbi_selection_outputs"
output_dir.mkdir(exist_ok=True)

df = pd.read_csv(input_file, dtype=str)
df.columns = [c.strip() for c in df.columns]

def find_col(possible_names, columns):
    lower_map = {c.lower(): c for c in columns}
    for name in possible_names:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None

organism_col = find_col([
    "Organism Name", "Organism Scientific Name", "organism_name",
    "Species", "species", "Organism"
], df.columns)

accession_col = find_col([
    "Assembly Accession", "assembly_accession", "Accession",
    "accession", "Assembly"
], df.columns)

assembly_level_col = find_col([
    "Assembly Level", "assembly_level", "Level", "level"
], df.columns)

scaffold_count_col = find_col([
    "Assembly Stats Number of Scaffolds", "Scaffold Count", "scaffold_count",
    "Scaffold count", "Scaffolds", "scaffolds", "Number of scaffolds"
], df.columns)

release_date_col = find_col([
    "Assembly Release Date", "Release Date", "release_date",
    "Submission Date", "submission_date", "Date", "date"
], df.columns)

gcf_col = find_col([
    "GCF_ID", "gcf_id", "RefSeq", "refseq_accession"
], df.columns)

gca_col = find_col([
    "GCA_ID", "gca_id", "GenBank", "genbank_accession"
], df.columns)

if organism_col is None:
    raise ValueError("Cannot detect organism/species name column.")

def normalize_species_name(name):
    if pd.isna(name):
        return None
    name = str(name).strip()
    name = re.sub(r"\s+", " ", name)

    parts = name.split()
    if len(parts) < 2:
        return name

    if len(parts) >= 2 and parts[1].lower() == "sp.":
        return name

    return " ".join(parts[:2])

df["source_organism_name"] = df[organism_col]
df["normalized_name"] = df[organism_col].apply(normalize_species_name)

def infer_accession(row):
    for col in [accession_col, gcf_col, gca_col]:
        if col is not None:
            val = row.get(col, None)
            if pd.notna(val) and str(val).strip() != "":
                return str(val).strip()
    return ""

df["best_accession"] = df.apply(infer_accession, axis=1)

def is_refseq(acc):
    acc = str(acc)
    return 1 if acc.startswith("GCF_") else 0

df["is_refseq"] = df["best_accession"].apply(is_refseq)

assembly_rank_map = {
    "complete genome": 1,
    "chromosome": 2,
    "scaffold": 3,
    "contig": 4
}

def assembly_rank(x):
    if pd.isna(x):
        return 999
    x = str(x).strip().lower()
    return assembly_rank_map.get(x, 999)

if assembly_level_col is not None:
    df["assembly_rank"] = df[assembly_level_col].apply(assembly_rank)
else:
    df["assembly_rank"] = 999

def to_number(x):
    try:
        return float(str(x).replace(",", ""))
    except:
        return 1e12

if scaffold_count_col is not None:
    df["scaffold_count_num"] = df[scaffold_count_col].apply(to_number)
else:
    df["scaffold_count_num"] = 1e12

if release_date_col is not None:
    df["release_date_parsed"] = pd.to_datetime(df[release_date_col], errors="coerce")
else:
    df["release_date_parsed"] = pd.NaT

df["release_date_sort"] = df["release_date_parsed"]

sort_cols = ["normalized_name", "is_refseq", "assembly_rank", "scaffold_count_num", "release_date_sort"]
ascending_flags = [True, False, True, True, False]

df_sorted = df.sort_values(by=sort_cols, ascending=ascending_flags).copy()
representatives = df_sorted.groupby("normalized_name", as_index=False).first()

rep_keys = set(
    representatives["normalized_name"].astype(str) + "||" +
    representatives["best_accession"].astype(str)
)

df_sorted["rep_key"] = df_sorted["normalized_name"].astype(str) + "||" + df_sorted["best_accession"].astype(str)
variants = df_sorted[~df_sorted["rep_key"].isin(rep_keys)].copy()

representatives["selection_rule"] = (
    "Selected as representative by priority: RefSeq > assembly level > scaffold count > release date"
)

variants["selection_rule"] = (
    "Not selected as representative; retained as intraspecies variant"
)

all_out = output_dir / "ncbi_cleaned_with_normalized_names.tsv"
rep_out = output_dir / "ncbi_representatives_candidate.tsv"
var_out = output_dir / "ncbi_intraspecies_variants_candidate.tsv"

df_sorted.to_csv(all_out, sep="\t", index=False)
representatives.to_csv(rep_out, sep="\t", index=False)
variants.to_csv(var_out, sep="\t", index=False)

print("Total input rows:", len(df))
print("Representative rows:", len(representatives))
print("Variant rows:", len(variants))
print("Check total:", len(representatives) + len(variants))
print(all_out)
print(rep_out)
print(var_out)
