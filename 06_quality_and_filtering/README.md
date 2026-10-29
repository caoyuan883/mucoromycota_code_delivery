# 06 Quality assessment and post-filtering

**Purpose.** BUSCO (genome and protein mode) and OMArk assessment; OMArk-based
contig audit and post-filtering of 10 genomes; comparison of 12 genomes with
their public annotations.

**Dependencies.** BUSCO 6.0.0 (offline lineages `fungi_odb12`,
`mucoromycota_odb12`), OMArk 0.5.0 / OMAmer (`LUCA.h5`), DIAMOND (reference
comparison), Python 3 with pandas where imported. See `docs/DEPENDENCIES.md`.

## A. BUSCO and OMArk tables

| script | configuration | notes |
|---|---|---|
| `run_busco_genome_offline.sh` | options in the script header | genome mode |
| `build_genome_qc_master_table.py` | `--help` | genome QC master table |
| `run_busco_protein_mucoromycota.sh` | options in the script header | protein mode |
| `run_busco_213_longest_isoform.sh` | env `PROTEIN_DIR`, `PROJECT`, `BUSCO_DOWNLOAD_DIR` | representative proteomes of the processing batch |
| `make_protein_qc_master_table.py` | `--protein-dir --busco-mucoromycota --busco-fungi --braker2-dir --braker3-dir --out-prefix` | |
| `run_omark_212.sh` | env `PROT_DIR`, `OMAMER_DB`, `OUT_DIR` | |
| `collect_omark_results.py` | `--root --out-prefix`; optional `--status`, `--expected 212` | |

Some historical scripts retain batch-specific sample counts. The released
dataset comprises the 202 genome records listed in Workbook 3.

## B. Contig audit and post-filtering (`omark_filtering/`)

Order: `audit_omark_contigs.py` -> `filter_omark_contigs.py` ->
`workbook7_sheet2_summarize_filtered_annotations.py`.

```bash
python3 omark_filtering/audit_omark_contigs.py --omark-dir OMARK --gtf-dir PREFILTER_GTF \
        --genome-dir GENOMES --genome-suffix .fna --out V5_OUT
python3 omark_filtering/filter_omark_contigs.py --v5-dirs V5_ACAU_OUT V5_OUT \
        --gtf-dir PREFILTER_GTF --faa-dirs PREFILTER_FAA --genome-dir GENOMES --out V6_OUT
```

`--faa-dirs` accepts several folders; `--allow-skip` lets `filter_omark_contigs.py`
finish when samples are skipped.

* Rule RU (`filter_omark_contigs.py`): among contigs carrying OMArk
  contaminant genes, a contig is removed if its GC content lies outside the
  host P5-P95 range and it either carries no consistent or inconsistent genes
  (FO), or carries consistent genes (FA) with a contaminant share >= 0.8 and a
  GC content closer to the FO median than to the host median.
* For *Acaulospora colombiana*, an additional sample-specific filter removed
  contigs with GC content >=43%, based on the separation of the high-GC
  component in the bimodal GC distribution and the associated OMArk non-target
  taxonomic signal. In the code this is `GC_EXCEPTION = {"Acaulospora_colombiana": 43.0}`
  (GC = 100 x (G+C)/(A+C+G+T), ambiguous bases excluded); it applies to all
  contigs of this assembly and to no other sample.
  The run for this sample used `historical/omark_filtered_v6/omark_contam_v5_acau.py`,
  which differs from `audit_omark_contigs.py` only in its input and output paths.
* `--v5-dirs` order sets priority (the *Acaulospora* run first).
* `--gtf-dir` must point to pre-filter GTF files, e.g. a view from
  `05_representative_isoforms/build_prefilter_view.py`.
* `filter_omark_contigs.py` falls back to `<sample>*.gtf` (prefix match) when
  the exact file name is absent.
* `audit_omark_contigs.py` refuses a non-empty output folder and writes
  `00_skipped_samples.tsv`; `filter_omark_contigs.py` exits 3 if samples were
  skipped (unless `--allow-skip`).
* Result of the released run: 10 genomes filtered, 20,927 proteins removed
  (`evidence/omark_filter_summary_10.tsv`); labels in `inputs/filtered10_labels.txt`.

## C. Public reference comparison (`reference12_validation/`)

Compares the annotation of 12 genomes with their public annotations. The 12
records, their source identifiers (with version), BRAKER system, subphylum and
public annotation source are listed in `inputs/reference12_public_annotations.tsv`.
Public annotation files are not included; they are obtained from the listed
sources. All scripts read
`REF12_BASE` (module folder with `12_faa_longest_isoform/`,
`12_faa_own_longest_isoform/`, `12_gff/`, `12_gtf_own/`).

Order: `make_longest_isoform_from_gff.py` -> `run_busco_ref12.sh`
(`REF12_BUSCO_PROJECT`, `BUSCO_DOWNLOAD_DIR`) -> `collect_busco_ref12.py` ->
`run_protein_comparison.sh` (optional argument: one strain) -> `compare_by_coordinates.py` ->
`gene_structure_stats.py`. `make_longest_isoform_from_gff.py` applies its own
longest-isoform rule to public GFF files; it is not the stage 05 rule.

**Outputs.** QC master tables, OMArk summary, filtered GTF/FASTA/FAA,
removed-gene lists, reference comparison tables.

**Historical copies.** `historical/` holds BUSCO runner variants, plotting and
spot-check scripts and `_versions/` (several files with the same name and
different content; their order is not recorded).

## Script names in run logs

The run logs and `historical/` copies use the earlier names below (full list:
`docs/RENAME_MAP.tsv`). Output folders, model names and data labels keep their
original names.

| release path | historical path | note |
|---|---|---|
| `build_genome_qc_master_table.py` | `04_build_genome_qc_master_table.py` |  |
| `omark_filtering/audit_omark_contigs.py` | `omark_filtered_v6/omark_contam_v5.py` |  |
| `omark_filtering/filter_omark_contigs.py` | `omark_filtered_v6/omark_filter_v6.py` |  |
| `reference12_validation/collect_busco_ref12.py` | `12_faa_12_gff_technical_validation/collect_busco_ref12.py` |  |
| `reference12_validation/compare_by_coordinates.py` | `12_faa_12_gff_technical_validation/compare_by_coordinates.py` |  |
| `reference12_validation/gene_structure_stats.py` | `12_faa_12_gff_technical_validation/gene_structure_stats.py` |  |
| `reference12_validation/make_longest_isoform_from_gff.py` | `12_faa_12_gff_technical_validation/make_longest_isoform_from_gff.py` |  |
| `reference12_validation/run_busco_ref12.sh` | `12_faa_12_gff_technical_validation/run_busco_ref12.sh` |  |
| `reference12_validation/run_protein_comparison.sh` | `12_faa_12_gff_technical_validation/run_protein_comparison.sh` |  |
| `run_busco_genome_offline.sh` | `03_run_busco_genome_offline.sh` |  |
| `run_busco_protein_mucoromycota.sh` | `04_run_busco_protein_mucoromycota.sh` |  |
