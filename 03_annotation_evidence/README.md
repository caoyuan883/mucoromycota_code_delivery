# 03 RNA-seq evidence for BRAKER

**Purpose.** Find public RNA-seq runs, select paired-end runs, download, QC,
align, merge per species and decide the BRAKER mode per genome.

**Inputs.** Species list, SRA run metadata, genome paths (user-supplied).

**Dependencies.** SRA Toolkit (esearch/efetch), fastp, HISAT2, SAMtools,
Python 3 with pandas/openpyxl (versions in `docs/DEPENDENCIES.md`).

**Configuration.** `check_sra.sh`, `get_srr_info.sh` and `filter_paired_runs.py`
read and write fixed file names in the current folder (`species_list.txt` ->
`sra_availability.csv` -> `sra_detailed_runinfo.csv` -> `sci_grade_download_list_v2.csv`);
run them in a dedicated folder. `map_genome_paths.py` and
`build_species_level_plan.py` take `--input/--output`.

```bash
bash download_qc_align.sh RUN_TABLE.xlsx WORK_DIR
bash summarize_alignment_qc.sh RNA_SEQ_DIR
bash merge_species_bam.sh RUN_BAM_DIR MERGED_BAM_DIR
bash map_species_braker_mode.sh GENOME_DIR MERGED_BAM_DIR MAP_OUT_CSV
```

`download_qc_align.sh` downloads the selected SRA runs (SRA Toolkit) into
`WORK_DIR`; both arguments have defaults in the script
(`18_tables_sra_run_download_list_with_genome_path.xlsx`, `rnaseq_braker_work`).

**Order.** `check_sra.sh` -> `get_srr_info.sh` -> `filter_paired_runs.py` -> `map_genome_paths.py` ->
`build_species_level_plan.py` -> `download_qc_align.sh` -> `summarize_alignment_qc.sh` ->
`merge_species_bam.sh` -> `map_species_braker_mode.sh`.

**Outputs.** Per-species merged BAM, alignment QC summary, `species_mode_mapping` CSV.

**Mode rule.** `map_species_braker_mode.sh` assigns BRAKER3 when a merged BAM
exists and BRAKER2 otherwise.

**Historical copies.** `historical/` holds the 12-JGI sub-batch variants and the
originals of the refactored scripts.

## Script names in run logs

The run logs and `historical/` copies use the earlier names below (full list:
`docs/RENAME_MAP.tsv`). Output folders, model names and data labels keep their
original names.

| release path | historical path | note |
|---|---|---|
| `build_species_level_plan.py` | `12_build_species_level_plan.py` |  |
| `check_sra.sh` | `08_check_sra.sh` |  |
| `download_qc_align.sh` | `13_download_qc_align_pipeline.sh` |  |
| `filter_paired_runs.py` | `10_filter_paired_runs.py` |  |
| `get_srr_info.sh` | `09_get_srr_info.sh` |  |
| `map_genome_paths.py` | `11_map_genome_paths.py` |  |
| `map_species_braker_mode.sh` | `16_map_species_braker_mode.sh` |  |
| `merge_species_bam.sh` | `15_merge_species_bam.sh` |  |
| `summarize_alignment_qc.sh` | `14_summarize_alignment_qc.sh` |  |
