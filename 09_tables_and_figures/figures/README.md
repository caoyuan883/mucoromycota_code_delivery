# 09 Figures - initial plotting scripts

These scripts generate the initial plots used for figure preparation.
The publication figures were subsequently refined manually in
Adobe Illustrator 2026. The scripts reproduce the initial plotting
stage rather than the exact final appearance of the publication figures.

## Files

| file | content |
|---|---|
| `figure_index.tsv` | one row per plotting script: topic, input data, required fields, outputs, manuscript figure, validation status |
| `plot_sources.tsv` | original file name and SHA256 -> repository path and SHA256, label in the original file name and in the script header, changes made |
| `historical/` | the original scripts, byte-identical, under their original file names (including the trailing `-N` of some names) |
| `historical/diffs/` | unified diff between each original script and its repository copy |

`manuscript_figure` is `UNASSIGNED` for every script. The labels in the original
file names (Fig1a, Fig1b, Fig2, FIG3, FIG4, Fig5, FIG7-12, FIG13, Fig14) are kept
in `plot_sources.tsv` only and are not figure numbers; some differ from the
script headers (for example `Fig14_...` has the header "Fig. 2").

## Script -> plot draft

| release script | plot draft | language / packages |
|---|---|---|
| `plot_database_overlap_annotation_availability.R` | NCBI/JGI overlap circles and public-protein availability bar | R: ggplot2 |
| `plot_busco_completeness_by_source.R` | BUSCO completeness of genome, public protein and BRAKER protein sets by subphylum | R: ggplot2, dplyr, tidyr, readr, readxl, patchwork |
| `plot_workflow_diagram.R` | workflow scheme | R: ggplot2, tibble; optional ragg, svglite |
| `plot_umbelopsis_ramanniana_busco.R` | BUSCO bars for *Umbelopsis ramanniana* | R: ggplot2, dplyr, tidyr, readr |
| `plot_genome_species_tree.R` | circular genome-based (G) tree | R: ape, ggplot2, dplyr, readr, stringr, Bioconductor treeio, ggtree |
| `plot_omark_contamination_filtering.py` | OMArk before/after filtering for the 10 filtered records (panels a, b) | Python 3: pandas, matplotlib |
| `plot_busco_genome_vs_protein_by_subphylum.R` | BUSCO genome vs protein, A4 pages per subphylum | R: ggplot2, dplyr, tidyr, readxl, readr, grid |
| `plot_busco_genome_reference_study_12.R` | BUSCO comparison for the 12 strains with public reference annotations | R: base graphics + grid; readxl only for Excel input |
| `plot_protein_tree_annotation.R` | protein-based (P) tree with annotation coverage bars and violin panel | R: ape, ggplot2, dplyr, readr, stringr, ggnewscale, patchwork, scales, Bioconductor treeio, ggtree |

Older plotting scripts are kept unchanged in their stage folders and listed in
`figure_index.tsv` as historical scripts:
`06_quality_and_filtering/historical/plot_*.r` and
`09_tables_and_figures/historical/Mucoromycota_tree_plot/plot_tree.py`.
`FIG6_build_omark_202-7.py` is byte-identical to
`06_quality_and_filtering/historical/build_omark_202.py`; it builds a table, not a
plot (see `plot_sources.tsv`).

## Inputs

| script | input | produced by |
|---|---|---|
| database overlap | none; the counts are constants in the script | - |
| BUSCO by source | per-record table with `subphylum, species_no, species_name_current, genome_BUSCO_C, public_protein_BUSCO_C, braker_protein_BUSCO_C` | stage 06 BUSCO results |
| workflow | none | - |
| *U. ramanniana* BUSCO | `Umbelopsis_ramanniana_BUSCO_completeness*.tsv` | stage 06 BUSCO results |
| G tree | `G_*.treefile` + `tip_metadata.tsv` | stage 08 tree; tip metadata table |
| OMArk filtering | `contamination_filtering_before_after_10strains_v2.tsv` | stage 06 OMArk before/after results |
| BUSCO G vs P | `longest_isoform_busco_genome_protein_mucoromycota_202_v3-2.xlsx`, sheet `1_简表_202株` | stage 06 longest-isoform BUSCO |
| 12-strain BUSCO | `busco_12strains_GTR_long_v1.tsv` | stage 06 `reference12_validation/` |
| P tree | folder with `P_*.treefile`, `Fig2_tip_data.tsv` (+ optional conflict tables) | stage 08 tree and comparison |

Required columns are listed in `figure_index.tsv`; the scripts stop when the
input or a required column is missing.

## Commands

Run from any folder; write into a new output folder.

```bash
Rscript 09_tables_and_figures/figures/plot_database_overlap_annotation_availability.R --outdir=OUT
Rscript 09_tables_and_figures/figures/plot_workflow_diagram.R --outdir=OUT
Rscript 09_tables_and_figures/figures/plot_busco_completeness_by_source.R --input=TABLE.tsv --outdir=OUT
Rscript 09_tables_and_figures/figures/plot_umbelopsis_ramanniana_busco.R --input=Umbelopsis_ramanniana_BUSCO_completeness.tsv --outdir=OUT
Rscript 09_tables_and_figures/figures/plot_genome_species_tree.R --tree=G_occ80_g90.treefile --metadata=tip_metadata.tsv --outdir=OUT
python3 09_tables_and_figures/figures/plot_omark_contamination_filtering.py --input=contamination_filtering_before_after_10strains_v2.tsv --outdir=OUT
Rscript 09_tables_and_figures/figures/plot_busco_genome_vs_protein_by_subphylum.R --input=longest_isoform_busco_genome_protein_mucoromycota_202_v3-2.xlsx --outdir=OUT
Rscript 09_tables_and_figures/figures/plot_busco_genome_reference_study_12.R --input=busco_12strains_GTR_long_v1.tsv --outdir=OUT
Rscript 09_tables_and_figures/figures/plot_protein_tree_annotation.R --input-dir=DIR_WITH_P_TREE_AND_TIP_TABLE --outdir=OUT
```

Each script prints its usage with `--help`. Without options the scripts behave as
before (RStudio file/folder chooser or the default paths in the script). The scripts
do not install packages; a missing package stops the script with its name.

## Output formats

PDF for all R scripts; additionally PNG (most), TIFF (database overlap, BUSCO by
source), SVG (workflow, *U. ramanniana*); the Python script writes PDF, SVG and PNG.
Several scripts also write small check tables (CSV/TSV). Fonts in PDF/PNG depend on
the fonts installed on the machine.

## Changes against the original scripts

Only the command-line options above, `--help`, and the package check (automatic
installation replaced by an error message) were changed; each change is marked
`[release]` in the code and shown in `historical/diffs/`. Data reading, data
transformations, filters, ordering, tree handling, support values, colours, layout
and styling are as in the original scripts.
