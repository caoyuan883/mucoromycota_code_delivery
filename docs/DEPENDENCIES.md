# Dependencies

Versions are listed with their source: "env" = `10_release_packaging/evidence/software_versions.txt`
(versions collected in the original environment); "log" = header of a tool
output or run log; "file" = a version file or path named in the column; "Workbook 4" = software table of the dataset.

| tool | stage | version | source | note |
|---|---|---|---|---|
| RepeatModeler | 02 | 2.0.7 | env | |
| RepeatMasker | 02 | 4.2.3 | env | |
| RMBlast | 02 | 2.14.1+ | env | |
| HISAT2 | 03 | 2.2.1 | env (GeneMark-ETP bundled copy) | |
| fastp | 03 | 0.23.4 | Workbook 4 | |
| SRA Toolkit | 03 | 3.4.1 | Workbook 4 | |
| SAMtools | 03/04 | 1.16.1 and 1.13 | env (two binaries) | |
| BRAKER (braker.pl) | 04 | 3.0.8 | env | one script tree for both batch systems (`$ENV_ROOT/braker2_bundle_transfer`) |
| AUGUSTUS | 04 | 3.5.0 | env | |
| ProtHint | 04 | 2.6.0 | env | |
| DIAMOND | 04 | 0.9.24 and 2.1.24 | env | `config.sh` uses `third_party/diamond/diamond_v2` |
| compleasm | 04 | 0.2.8 | env | |
| TSEBRA `best_by_compleasm.py` | 04 | 1.0.1 with local modification | see below | the modified file is required |
| gffread | 04/05 | 0.12.7 | env | |
| gtf2gff.pl (AUGUSTUS scripts) | 04/05 | revision not recorded | - | |
| BUSCO | 06 | 6.0.0 | Workbook 4 | lineages `fungi_odb12`, `mucoromycota_odb12` |
| OMArk | 06 | 0.5.0 | Workbook 4 | OMAmer database `LUCA.h5` |
| InterProScan | 07 | 5.78-109.0 (InterPro release 109) | file: `07_functional_annotation/evidence/interproscan_VERSIONS.txt`; installation path in the run scripts | |
| InterProScan member databases | 07 | Pfam 38.2, PANTHER 19.0, SUPERFAMILY 1.75, CDD 3.21, SMART 9.0, PROSITE 2026_01, AntiFam 8.0 | same file | |
| eggNOG-mapper | 07 | 2.1.13 | log (annotation file headers) | |
| MAFFT | 08 | 7.526 | env; header of `align_trim_concat.py` | |
| ClipKIT | 08 | 2.13.2 | header of `align_trim_concat.py` | |
| IQ-TREE | 08 | 3.1.0 | log: `08_phylogenomics/evidence/iqtree_logs/` (all four trees) | per-tree settings in `config/phylogenomics_logged_runs.tsv` |
| DendroPy | 08 | 5.1.0 | tests | only for `compare_four_trees.py` |
| R | 02/06/09 | 4.4.1 | env | historical plotting scripts: ggplot2, dplyr, readr, tidyr |
| R packages for `09_tables_and_figures/figures/` | 09 | versions not recorded | - | ggplot2, dplyr, tidyr, readr, readxl, patchwork, tibble, stringr, scales, ape, ggnewscale; Bioconductor treeio, ggtree; optional ragg, svglite, rstudioapi |
| Python plotting (`plot_omark_contamination_filtering.py`) | 09 | versions not recorded | - | pandas, matplotlib |
| Python | all | 3.10 / 3.13 | env | repository scripts tested with Python 3.14 |
| ete3 | 01 | - | imports | `NCBITaxa` with a local taxonomy database |
| pandas, openpyxl, python-docx | 01/03/06/09 | - | imports | only for the scripts that import them |

## best_by_compleasm.py

The BRAKER gene sets were post-processed with `best_by_compleasm.py` from the
TSEBRA copy in the BRAKER environment (TSEBRA 1.0.1). This study used a locally
modified version of the script, not the standard upstream file; the
modification skips HMM hits without a score cutoff.

| revision | SHA256 | note |
|---|---|---|
| production (skip HMM hits without a cutoff) | `b51574065bc255e424a2eadf713b7905dbbe2e1c168177293e9c99f929c99f4d` | used for the runs reported here; server backup name `best_by_compleasm.py.bak_1784124955` |
| later (identifier-mapping fix) | `e01e08daf8fc052747b1319fcdc73a84f321b54272914ecefa20930f1c1e1cfd` | listed in `10_release_packaging/evidence/software_versions.txt`; not the production revision |

The modified file is required to run this step. It is not included in this
repository and is not provided by standard BRAKER or TSEBRA installations. A
checksum identifies a revision but does not replace its source. Method details:
`04_structural_annotation/README.md`.

## Databases and test requirements

Large reference databases (BUSCO lineages, OMAmer `LUCA.h5`, InterProScan data,
eggNOG database, BRAKER protein evidence) must be installed before the
corresponding stages are run. The RNA-seq preparation scripts of stage 03 can
download the selected SRA runs as documented in `03_annotation_evidence/README.md`.

The default test run needs only the Python standard library. DendroPy
(`pip install dendropy==5.1.0`) is needed for the tree-comparison tests; without
it these tests are reported as NOT_RUN. See `tests/README.md`.
