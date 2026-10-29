# Small inputs

Small tables used by the release scripts and the tests. Do not edit them in
place; if a value has to change, add a new version and record the reason in
`docs/CHANGES.md`.

| file | content | origin |
|---|---|---|
| `Workbook3_Final_202.tsv` | TSV export of Workbook 3 (202 records: Record No., source identifier, names, subphylum) | curated table, exported unchanged |
| `Workbook3_Naming_notes_6.tsv` | notes for the six records that keep the source name | exported unchanged |
| `table7_species_name_mapping_202.txt` | name mapping used for Table 7 | exported unchanged |
| `02_per_strain_functional_summary_final.tsv` | per-genome functional statistics of the release (202 rows, column `filtered`) | output of `07_functional_annotation/recompute_filtered10_functional_stats.py` |
| `per_strain_release.tsv` | per-genome counts of the independent recheck | output of `07_functional_annotation/recheck_func_numbers.py` |
| `release_labels_202.txt` | the 202 release labels, one per line | `strain` column of `per_strain_release.tsv` |
| `filtered10_labels.txt` | the 10 post-filtered labels | rows with `filtered = yes` in `02_per_strain_functional_summary_final.tsv` |
| `label_overrides_6.tsv` | release labels of records 41, 94, 111, 127, 151, 152 (source name retained) | naming exception (`docs/DATA_LINKAGE.md`); optional cross-check with `derive_subphylum_map.py --overrides` |
| `reference12_public_annotations.tsv` | the 12 genomes compared with public annotations: Record No., release label, source identifier with version, BRAKER system, subphylum, public annotation source | Workbook 3 and the public annotation records |

## Workbook 3 export

| item | value |
|---|---|
| workbook file | `03_final_dataset-2.xlsx` |
| SHA256 | `d0595481587f47330057725262dcd4a6b6e6ef435cd6db84e2345d438cf72607` |
| `Workbook3_Final_202.tsv` | sheet `Final_202`, header row 4, data rows 5-206 (202 records) |
| `Workbook3_Naming_notes_6.tsv` | sheet `Naming_notes_6`, header row 4, data rows 5-10 |

Both sheets are exported cell by cell without changes (title and note rows 1-3
are omitted). Workbook 3 is the reference for record selection and names. In
this export, the column "Name used in this study" of records 41, 94, 111, 127,
151 and 152 gives the source name, which is also their release label.

See `docs/DATA_LINKAGE.md` for how these tables connect to the per-genome files.
