# Usage notes

## Large inputs

Genome assemblies, RNA-seq alignments, BUSCO lineages, the OMAmer database,
InterProScan data, the eggNOG database and the BRAKER protein database are not
part of this repository. The large reference databases must be installed
before running the corresponding stages; their locations are set in the
configuration (`config/*.example`, environment variables or command-line
options). The RNA-seq preparation scripts of stage 03 query NCBI SRA and
download the selected runs as described in `03_annotation_evidence/README.md`.

## Identifiers

Gene and transcript IDs (`g1`, `g1.t1`) are unique only within one genome.
Combine them with the release label before merging tables across genomes.
Genome records are identified by `Record No.` and source identifier
(Workbook 3); never merge genomes by their current species name.

## Pre-filter and release data

Ten genomes were post-filtered after the OMArk contig audit. The release
contains 2,618,164 representative proteins; the pre-filter state contains
2,639,091. Statistics that refer to the pre-filter state must be computed on a
pre-filter view (`05_representative_isoforms/build_prefilter_view.py`), and the
two states must not be mixed in one table.

## Sample counts in historical scripts

Some historical scripts retain batch-specific sample counts. The released
dataset comprises the 202 genome records listed in Workbook 3.

## Trees

Comparing the four trees needs only their tree files
(`08_phylogenomics/compare_four_trees.py`). Tree reconstruction requires the
corresponding supermatrices and locus lists; these are not bundled with this
repository. The IQ-TREE settings of the four trees, taken from their run logs
(`08_phylogenomics/evidence/iqtree_logs/`), are in
`config/phylogenomics_logged_runs.tsv`.

## Figures

The plotting scripts in `09_tables_and_figures/figures/` produce initial plots.
The publication figures were refined manually in Adobe Illustrator 2026
(see `09_tables_and_figures/figures/README.md`).

## Historical copies

Files under `historical/` keep the original workstation paths and names. Use the
top-level scripts of each stage; `docs/PROVENANCE.tsv` and `docs/RENAME_MAP.tsv`
link them to the historical copies and run logs.
