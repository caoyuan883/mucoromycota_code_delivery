#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script : annotate_ncbi_taxonomy.py
Purpose: Query the NCBI Taxonomy database using organism names to extract 
         hierarchical taxonomic ranks (subphylum to species).
"""

import argparse
import sys
import re
import pandas as pd
from pathlib import Path
try:
    from ete3 import NCBITaxa
except ImportError:
    sys.exit("ERROR: The 'ete3' package is required. Please install it (e.g., pip install ete3).")

def normalize_name(name):
    """
    Cleans the input name by removing underscores, extra spaces, 
    and standardizing the 'sp.' formatting.
    """
    s = str(name).strip()
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"\bsp\.\s*", "sp. ", s)
    s = re.sub(r"\bsp\s+", "sp. ", s)
    s = s.strip()
    return s

def candidate_names(name):
    """
    Generates candidate names for NCBI query. If 'sp.' is present, 
    it also attempts to query just the genus as a fallback.
    """
    s = normalize_name(name)
    cands = [s]

    if " sp. " in s:
        genus = s.split(" sp. ")[0].strip()
        cands.append(genus)

    # Return unique candidates preserving order
    return list(dict.fromkeys(cands))

def query_taxonomy(ncbi, name):
    """
    Queries NCBI for the given name candidates and extracts taxonomy ranks.
    """
    for cand in candidate_names(name):
        try:
            trans = ncbi.get_name_translator([cand])
            if cand in trans and trans[cand]:
                taxid = trans[cand][0]
                lineage = ncbi.get_lineage(taxid)
                ranks = ncbi.get_rank(lineage)
                names = ncbi.get_taxid_translator(lineage)

                rank2name = {}
                for tid in lineage:
                    r = ranks.get(tid, "")
                    n = names.get(tid, "")
                    if r and n:
                        rank2name[r] = n

                genus = rank2name.get("genus", "")
                species = rank2name.get("species", "")

                # Fallback inferences if NCBI lacks strict rank assignments
                if not genus:
                    parts = cand.split()
                    if len(parts) >= 1:
                        genus = parts[0]

                if " sp. " in cand:
                    species = cand
                elif not species:
                    parts = cand.split()
                    if len(parts) >= 2:
                        species = " ".join(parts[:2])

                return {
                    "subphylum": rank2name.get("subphylum", ""),
                    "class": rank2name.get("class", ""),
                    "order": rank2name.get("order", ""),
                    "family": rank2name.get("family", ""),
                    "genus": genus,
                    "species": species
                }
        except Exception:
            pass

    # Complete fallback if no NCBI match is found
    s = normalize_name(name)
    parts = s.split()
    genus = parts[0] if len(parts) >= 1 else ""
    species = s if len(parts) >= 2 else ""
    return {
        "subphylum": "",
        "class": "",
        "order": "",
        "family": "",
        "genus": genus,
        "species": species
    }

def main():
    parser = argparse.ArgumentParser(description="Extract taxonomy ranks for organisms using ETE3 NCBITaxa.")
    parser.add_argument("--input-csv", required=True, help="Input CSV file containing sample/organism names in the first column.")
    parser.add_argument("--out-csv", required=True, help="Output CSV file for the taxonomy results.")
    args = parser.parse_args()

    input_file = Path(args.input_csv).resolve()
    output_file = Path(args.out_csv).resolve()

    if not input_file.is_file():
        sys.exit(f"ERROR: Input file not found: {input_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Load input data
    df = pd.read_csv(input_file)
    if df.empty:
        sys.exit("ERROR: Input CSV is empty.")
        
    first_col = df.columns[0]
    samples = df[first_col].astype(str).tolist()

    # Initialize local NCBI Taxonomy database
    print("Initializing NCBI Taxonomy database (may take a moment to update/download)...")
    ncbi = NCBITaxa()

    print(f"Processing {len(samples)} samples...")
    rows = []
    for sample in samples:
        tax = query_taxonomy(ncbi, sample)
        rows.append({
            "sample": sample,
            "subphylum": tax["subphylum"],
            "class": tax["class"],
            "order": tax["order"],
            "family": tax["family"],
            "genus": tax["genus"],
            "species": tax["species"]
        })

    # Save output
    out = pd.DataFrame(rows, columns=["sample", "subphylum", "class", "order", "family", "genus", "species"])
    out.to_csv(output_file, index=False)

    print(f"\n===== Taxonomy Extraction Complete =====")
    print(f"Successfully generated: {output_file}")
    print("\nPreview of first 5 rows:")
    print(out.head(5).to_string(index=False))

if __name__ == "__main__":
    main()
