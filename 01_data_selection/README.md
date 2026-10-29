# 01 Data selection

**Purpose.** Select one representative genome assembly per taxon from NCBI and
JGI MycoCosm metadata, normalise names and annotate NCBI taxonomy. The curated
result is Workbook 3, the list of the 202 released genome records.

**Inputs.** NCBI datasets metadata (TSV) and the JGI MycoCosm table (CSV),
downloaded by the user.

**Dependencies.** Python 3, pandas; ete3 (`NCBITaxa` with a local taxonomy
database) for `annotate_ncbi_taxonomy.py`.

**Configuration.** Command-line options (`--help`).

**Order.**
1. `select_ncbi_representative_genomes.py`
2. `count_jgi_species.py`
3. `annotate_ncbi_taxonomy.py`
4. Manual curation into Workbook 3 (export: `inputs/Workbook3_Final_202.tsv`).

**Outputs.** Candidate genome lists and taxonomy tables used to compile
Workbooks 1-3.

**Checks.** Workbook 3 contains 202 records: Mucoromycotina 123,
Mortierellomycotina 57, Glomeromycotina 22 (checked by
`10_release_packaging/derive_subphylum_map.py`). Later tables join to Workbook 3
by `Record No.` and source identifier.

**Historical copies.** `historical/` holds earlier versions (`_versions/`), the
naming audit (`Mucoromycota_tree_plot/`) and `07_build_species_list.py`, which
has a `.py` suffix but contains shell commands (kept unchanged; test S03 lists it).

## Script names in run logs

The run logs and `historical/` copies use the earlier names below (full list:
`docs/RENAME_MAP.tsv`). Output folders, model names and data labels keep their
original names.

| release path | historical path | note |
|---|---|---|
| `annotate_ncbi_taxonomy.py` | `05_annotate_ncbi_taxonomy.py` |  |
| `count_jgi_species.py` | `02_count_jgi_species.py` |  |
| `select_ncbi_representative_genomes.py` | `01_select_ncbi_representatives_genomes.py` |  |
