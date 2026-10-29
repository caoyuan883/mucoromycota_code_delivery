# mucoromycota-annotation-scripts

Scripts used to build a structural and functional annotation resource for
202 Mucoromycota genome records (Mucoromycotina 123, Mortierellomycotina 57,
Glomeromycotina 22): genome selection, repeat masking, BRAKER2/BRAKER3 gene
prediction, representative isoforms, quality assessment and post-filtering,
functional annotation, phylogenomics and per-genome packaging.

## Data scope

* 202 genome records, listed with their source identifiers in Workbook 3.
* One representative transcript per gene: 2,618,164 representative proteins.
  Ten genomes were post-filtered after an OMArk contig audit (20,927 proteins
  removed; labels in `inputs/filtered10_labels.txt`).
* Five files per genome, packed as `<Subphylum>/<label>.tar.gz`:
  `<label>_longest_isoform.gtf`, `<label>_longest_isoform.gff3` (structural,
  converted from the GTF), `<label>_longest_isoform.faa`,
  `<label>_interproscan.tsv`, `<label>_emapper.annotations.tsv`.
* Four species trees from BUSCO single-copy loci (`fungi_odb12`):
  `P_occ80_g90`, `P_occ90_g90`, `G_occ80_g90`, `G_occ90_g90`
  (P = BUSCO on the annotated proteomes, G = BUSCO genome mode; occ = per-locus
  taxon occupancy; g90 = maximum column gap 0.9).

## Repository structure

```
01_data_selection/          genome selection (NCBI, JGI MycoCosm), taxonomy
02_repeat_masking/          RepeatModeler2 + RepeatMasker soft-masking, repeat tables
03_annotation_evidence/     public RNA-seq: selection, QC, alignment, BRAKER mode per genome
04_structural_annotation/   BRAKER2 (braker2/) and BRAKER3 (braker3/), two independent batch systems
05_representative_isoforms/ one representative transcript per gene; GTF -> GFF3
06_quality_and_filtering/   BUSCO, OMArk, contig audit and post-filtering, 12-genome reference comparison
07_functional_annotation/   InterProScan / eggNOG-mapper results, per-genome statistics
08_phylogenomics/           marker selection, matrices, IQ-TREE, four-tree comparison
09_tables_and_figures/      plotting scripts (figures/), table builders (historical/)
10_release_packaging/       release assembly, checks, per-subphylum packaging
config/                     example configuration files
inputs/                     small tables used by the scripts and tests
docs/                       dependencies, data linkage, usage notes, provenance, changes
tests/                      light-weight tests (synthetic data and small inputs)
```

In each stage folder, top-level scripts are the entry points. `historical/`
holds byte-identical earlier versions and the originals of refactored scripts
(`historical/pre_refactor/`); these copies keep the original workstation paths
and are kept for provenance. `evidence/` holds small logs and reports from the
original runs. The status of every file (`historical`, `refactored`,
`reconstructed`, `new_validation`, inputs and documentation), its source path and
SHA256 are listed in `docs/PROVENANCE.tsv`.

## Installation and configuration

The scripts need Python 3, bash and, per stage, the tools listed in
`docs/DEPENDENCIES.md`. Genome data and large reference databases are not
distributed; the databases must be installed beforehand and their locations
set by the user. The RNA-seq preparation scripts (stage 03) can download the
selected SRA runs as documented in `03_annotation_evidence/README.md`.

* Python and R entry points take paths as command-line options (`--help`).
* Shell entry points read environment variables and stop with a message such as
  `set REP_ROOT (...)` when one is missing.

```bash
cp config/paths.example.env config/paths.env              # edit paths
cp config/braker_site.env.example config/braker_site.env  # BRAKER site settings
set -a; source config/paths.env; set +a
export BRAKER_SITE_ENV=$PWD/config/braker_site.env
```

Scripts that create an output folder refuse to overwrite an existing one unless
`--overwrite` / `OVERWRITE=1` is given. `STRICT=1` in the release scripts turns
the expected totals (202 genomes, 2,618,164 proteins, 123/57/22) into hard checks.

## Main entry points

| stage | entry points |
|---|---|
| 01 | `select_ncbi_representative_genomes.py`, `count_jgi_species.py`, `annotate_ncbi_taxonomy.py` |
| 02 | `run_repeat_softmasking.sh`, `summarize_repeatmodeler_logs.sh`, `check_softmask.py` |
| 03 | `check_sra.sh` ... `map_species_braker_mode.sh` (order in the stage README) |
| 04 | `braker2/` and `braker3/`: `config.sh`, `preflight.sh`, `run_batch.sh`; `extract_braker_gtf_from_archives.py` |
| 05 | `filter_longest_isoforms_from_gtf.py`, `repair_3_longest_isoform_samples.py`, `gtf_to_gff3.sh`, `build_prefilter_view.py` |
| 06 | BUSCO/OMArk runners and tables; `omark_filtering/audit_omark_contigs.py`, `omark_filtering/filter_omark_contigs.py`; `reference12_validation/` |
| 07 | `ips_fix_20261011/` freeze chain; `summarize_functional_annotation.py`, `recompute_filtered10_functional_stats.py`, `assemble_functional_tables.sh`, `recheck_func_numbers.py` |
| 08 | `collect_busco_markers.py`, `diagnose_and_select_markers.py`, `align_trim_concat.py`, `matrix_qc_and_column_filter.py`, `run_iqtree_occ80.sh`, `compare_four_trees.py` |
| 09 | `figures/plot_*.R`, `figures/plot_omark_contamination_filtering.py` |
| 10 | `derive_subphylum_map.py`, `check_release_inputs.py`, `make_release_202.sh`, `pack_release_202.sh` |

Each stage `README.md` gives inputs, configuration, commands, outputs and checks.

### Examples using repository files only

```bash
# release label -> subphylum table used by pack_release_202.sh
python3 10_release_packaging/derive_subphylum_map.py \
    --workbook3 inputs/Workbook3_Final_202.tsv --labels inputs/release_labels_202.txt \
    --out /tmp/00_subphylum_map.tsv --record-map /tmp/record_label_map.tsv

# light-weight tests (minimal run)
python3 tests/run_tests.py --outdir /tmp/mucoro_tests
```

### Examples with external inputs

```bash
# compare the four trees (needs DendroPy and the four .treefile files)
python3 08_phylogenomics/compare_four_trees.py --tree-dir /path/to/trees --outdir /tmp/tree_comparison

# full test run with tree files and optional tools
python3 tests/run_tests.py --outdir /tmp/mucoro_tests --trees /path/to/trees \
    --gtf2gff /path/to/gtf2gff.pl --rscript /path/to/Rscript --plot-python /path/to/python3
```

## Data and workbook linkage

Genome records are identified by `Record No.` and source identifier in
Workbook 3; files use a release label (for example `Acaulospora_colombiana`).
Gene and transcript IDs (for example `g1.t1`) are unique only within one genome
and must be combined with the label. The mapping between records, labels,
subphyla and files is described in `docs/DATA_LINKAGE.md`.

## Usage notes

Short notes on large inputs, identifiers, pre-filter and release data, trees and
figures are in `docs/USAGE_NOTES.md`.

## Documentation

* `docs/DEPENDENCIES.md` - software, versions and the record each version is taken from
* `docs/DATA_LINKAGE.md` - records, labels, workbooks and per-genome files
* `docs/USAGE_NOTES.md` - practical notes for using the scripts and data
* `docs/PROVENANCE.tsv` - every file with source path, SHA256, status and change
* `docs/CHANGES.md` - changes made to the scripts for this repository
* `docs/RENAME_MAP.tsv` - current script names and the names used in the run logs
* `tests/README.md` - test groups and the reference test result

## License

MIT License; see `LICENSE`.
