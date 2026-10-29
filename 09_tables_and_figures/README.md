# 09 Tables and figures

**Layout.**
* `figures/` - plotting scripts for the initial figure drafts, `figure_index.tsv`
  and `plot_sources.tsv`; see `figures/README.md`.
* `historical/` - table builders kept with their original names: Table 6
  (`build_表6_*.py`), Table 7 (`build_table7.py`, `build_表7_*.py`), the
  genome-vs-protein BUSCO comparison
  (`make_longest_isoform_genome_protein_comparison.py`), the benchmark
  selection reports (`make_benchmark_*.py`, earlier copies in `_versions/`) and
  a tree plot (`Mucoromycota_tree_plot/plot_tree.py`). Several builders write
  to the same output path and their order is not recorded; run them with
  `--help` or read their headers before use.

**Workbooks.** The multi-sheet XLSX workbooks were assembled by copying table
values manually. Workbook 3 is the curated list of the 202 genome records
(export: `inputs/Workbook3_Final_202.tsv`); the links between the workbooks
and the scripts are listed in `docs/DATA_LINKAGE.md`.

**Dependencies.** Python 3 with pandas, openpyxl and python-docx for the table
builders; plotting packages are listed in `figures/README.md`.
