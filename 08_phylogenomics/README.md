# 08 Phylogenomics (four species trees)

**Purpose.** Build and compare four species trees from BUSCO single-copy loci
(lineage `fungi_odb12`): `P_*` from BUSCO on the annotated proteomes, `G_*`
from BUSCO genome mode; `occ80`/`occ90` = per-locus taxon occupancy of 0.80 or
0.90; `g90` = column gap filter 0.9. 207 taxa = 202 ingroup + 5 outgroups
(Neurospora_crassa, Aspergillus_nidulans, Coprinopsis_cinerea,
Conidiobolus_coronatus, Coemansia_reversa).

**Parameters.** `config/phylogenomics_logged_runs.tsv` lists, for each tree,
the settings recorded in its IQ-TREE log (`evidence/iqtree_logs/<tree>.log`):
IQ-TREE 3.1.0, input matrix `supermatrix_<tree>.fasta`, `-st AA`, model
JTT+F+G4, `-B 1000 -alrt 1000 -bnni`, seed 12345; occ80 trees `-T 20 -mem 80G`,
occ90 trees `-T 16 -mem 45G`. The executable path in the logs belongs to the
original workstation; the commands below take it from `IQ`.

**Dependencies.** Python 3; MAFFT 7.526 and ClipKIT 2.13.2 (or trimAl) for
`align_trim_concat.py`; IQ-TREE 3; DendroPy (tested with 5.1.0) for the comparison.

**Order.**
1. `collect_busco_markers.py --p-dir P_BUSCO_SC --g-dir G_BUSCO_SC --out ROOT`
   (QC, loci valid in both sources, occupancy tables).
2. `diagnose_and_select_markers.py --out ROOT` (locus lists per threshold).
3. `align_trim_concat.py --root ROOT --locus-list ROOT/reports/locus_list_occ80.txt --tag occ80` (optional `--trimmer clipkit`).
4. `matrix_qc_and_column_filter.py --root ROOT --tag occ80 --max-col-gap 0.9`.
5. `IQ=... WORK=... bash compare_substitution_models.sh`, then
   `bash run_iqtree_occ80.sh memcheck`, `bash run_iqtree_occ80.sh mf` and
   `bash run_iqtree_occ80.sh tree JTT+F+G4` (model from the `mf` stage); the
   `tree` stage runs `-B 1000 -alrt 1000 -bnni` as in the logged runs.
6. `compare_four_trees.py --tree-dir TREES --outdir OUT` (optional `--expected-ingroup 202`).

For occ90, run steps 3-4 with `--tag occ90` (locus list
`ROOT/reports/locus_list_occ90.txt`) and the tree command recorded in
the occ90 logs, once per side (`P`, `G`):

```bash
"$IQ" -s "$WORK/supermatrix_P_occ90_g90.fasta" -st AA -m JTT+F+G4 \
      -B 1000 -alrt 1000 -bnni -T 16 -mem 45G -seed 12345 --prefix P_occ90_g90
```

Tree reconstruction requires the corresponding supermatrices and locus lists;
these are not bundled with this repository. Comparing the four trees
(step 6) needs only the four `.treefile` files.

**Comparison output.** `pairwise_summary.tsv` and `differing_splits.tsv`.
Normalised RF = RF / (2n - 6) on the 202 ingroup tips after pruning the
outgroups; strong support = SH-aLRT >= 80 and UFBoot >= 95. On the released
tree files (test T01): P_occ80 vs G_occ80 share 192 of 199 splits, nRF 0.0352.

## Script names in run logs

The run logs and `historical/` copies use the earlier names below (full list:
`docs/RENAME_MAP.tsv`). Output folders, model names and data labels keep their
original names.

| release path | historical path | note |
|---|---|---|
| `align_trim_concat.py` | `step02_align_trim_concat_v2.py` |  |
| `collect_busco_markers.py` | `step01_collect_busco.py` |  |
| `compare_four_trees.py` | `compare_four_trees_dendropy.py` |  |
| `compare_substitution_models.sh` | `step03b_model_compare.sh` |  |
| `diagnose_and_select_markers.py` | `step01b_diagnose_and_select.py` |  |
| `matrix_qc_and_column_filter.py` | `step02c_matrix_qc_and_colfilter.py` |  |
| `run_iqtree_occ80.sh` | `step03_iqtree_occ80.sh` |  |
