# Changes against the original scripts

The scripts in stages 01-10 were adapted for use outside the original
workstations. The scientific rules (representative-isoform key and
tie-breaks, contamination rule and the *Acaulospora colombiana* GC filter,
protein cleaning, denominators, GO split, union/intersection, BUSCO lineages,
InterProScan member databases, matrix filters and substitution models) are
unchanged; the IQ-TREE support options of the public tree entry point were
aligned with the run logs (section 7). Per-file details with the SHA256 of the original and
the repository version: `PROVENANCE.tsv`. The original of every modified script
is kept under `<stage>/historical/pre_refactor/`; changed lines are marked
`[release refactor]` (or `[release]` in the plotting scripts).

## 1. Paths

Workstation paths (`/media/desk16/...`, `/home/py/Desktop/...`,
`/Users/caoyuan/Desktop/...`, `~/Desktop/...`) were replaced by command-line
options (Python, argparse) or required environment variables (shell,
`${VAR:?message}`). The original value is kept as a comment next to each
change. Implicit directory choices (first or last glob match) became explicit
arguments (`summarize_functional_annotation.py`, `recompute_filtered10_functional_stats.py`).

## 2. Input checks, overwrite protection, exit codes

* Missing input folders or files stop the script with a message and a non-zero exit.
* Scripts that write a new output folder refuse a non-empty one
  (`--overwrite` / `OVERWRITE=1` where re-packing is intended).
* `summarize_functional_annotation.py`, `recompute_filtered10_functional_stats.py`:
  non-zero exit when the run fails (the originals printed the exception and
  ended normally); `recompute_filtered10_functional_stats.py` exits 2 when a
  filtered file matches no table row.
* `filter_omark_contigs.py` exits non-zero when samples are skipped for missing
  GTF/FASTA (`--allow-skip` restores the original behaviour);
  `audit_omark_contigs.py` writes skipped samples to `00_skipped_samples.tsv`.
* Optional `STRICT=1` in `assemble_functional_tables.sh`,
  `subset_filtered10_functional_tables.sh`, `make_release_202.sh` and
  `pack_release_202.sh` turns the printed expectations into checks (exit 3).
* `set -euo pipefail` was not added globally; batch scripts keep their
  original error handling.

## 3. Portability

`gzcat` and `md5 -r` (macOS) were replaced by helpers that use `gzip -dc` and
`md5sum` when available (`md5 -r` otherwise). `md5sum` separates hash and name
with two spaces, `md5 -r` with one.

## 4. Checksum lists

`make_release_202.sh` and `pack_release_202.sh` both wrote `MD5SUMS.txt` into
the release folder, so packing overwrote the per-file checksums. The per-file
list is now `MD5SUMS_files.txt`; `MD5SUMS.txt` holds the archive checksums.

## 5. Script names and layout

* Scripts in stages 01-10 have functional snake_case names; stage numbers are
  kept. Map of old and new names: `RENAME_MAP.tsv`. Only references to the
  renamed scripts changed (`source`, `bash`, `python3` calls, `pgrep -f`
  patterns, usage lines, README/config/test paths). Run folders, deployed
  paths, AUGUSTUS species names, data labels, output prefixes and report names
  keep their original names.
* Folders: `04_structural_annotation/braker2/` and `braker3/`;
  `omark_filtered_v6/` -> `omark_filtering/`;
  `12_faa_12_gff_technical_validation/` -> `reference12_validation/`;
  the historical subfolder `ips_fix_20261011/fix_20261011_part2/` merged into `ips_fix_20261011/`.
* The BRAKER2 best_by_compleasm scripts (`06_bbc*.sh`, `09_rebuild_manifest.sh`)
  and `historical/10_rerun_bbc_serial.sh` keep their original names.
* `collect_busco_markers.py` moved from `06_quality_and_filtering/` to
  `08_phylogenomics/`, next to the other steps of the tree chain; content unchanged.
* `07_package_core.sh` (superseded by `package_results.sh`) moved unchanged to
  `04_structural_annotation/historical/braker2_batch_v3/`.

## 6. Plotting scripts

Nine plotting scripts were added to `09_tables_and_figures/figures/` with
`--help`, input/output options and an error message instead of automatic
package installation; originals and diffs are in `figures/historical/`. Data
handling and styling are unchanged.

## 7. Parameters aligned with the run logs

* `08_phylogenomics/run_iqtree_occ80.sh`: the tree stage adds `-bnni` and runs
  `-B 1000 -alrt 1000 -bnni`, as recorded in the four IQ-TREE run logs
  (`08_phylogenomics/evidence/iqtree_logs/`). The original script, without
  `-bnni`, is kept as `historical/pre_refactor/step03_iqtree_occ80.sh`.
* `04_structural_annotation/braker{2,3}/config.sh`: `THREADS_PER_JOB` (default 8)
  and `MAX_PARALLEL` (default 4) can be set from the environment
  (originals: 8 and 1 for BRAKER2, 8 and 4 for BRAKER3).
* `config/phylogenomics_logged_runs.tsv`: per-tree IQ-TREE settings taken from
  the four run logs (version 3.1.0, input matrix, model, support options,
  threads, memory, seed, run times), plus the occupancy (0.80/0.90) and
  column-gap (0.9) thresholds encoded in the matrix names. Threads and memory
  differ between the occ80 trees (20, 80G) and the occ90 trees (16, 45G).

## 8. Inputs and outputs

* `inputs/Workbook3_Final_202.tsv` was re-exported from the current Workbook 3
  (`03_final_dataset-2.xlsx`; source details in `inputs/README.md`). It differs
  from the earlier export only in the column "Name used in this study" of
  records 41, 94, 111, 127, 151 and 152, which now gives the source name. With
  this export `derive_subphylum_map.py` maps all 202 records without
  `--overrides`.
* `inputs/label_overrides_6.tsv`: labels of the same six records; kept as an
  explicit record and optional cross-check (`derive_subphylum_map.py --overrides`).
* `inputs/reference12_public_annotations.tsv`: the 12 genomes of the public
  reference comparison with source identifiers and public annotation sources.
* `07_functional_annotation/recheck_func_numbers.py` additionally writes
  `subphylum_medians.tsv` (unrounded and one-decimal medians). One `EXPECT`
  value changed: `median_any_Glo` 72.2 -> 72.3, the manuscript display value;
  the median is still recomputed (72.252888) and compared at one decimal.
  `evidence/recheck_report.txt` was written with the earlier value and is kept
  unchanged.
* `08_phylogenomics/evidence/iqtree_logs/`: the four IQ-TREE run logs.
* `LICENSE`: MIT License.

## 9. Reconstructed and new files

| file | status | purpose |
|---|---|---|
| `05_representative_isoforms/gtf_to_gff3.sh` | reconstructed | wrapper around the recorded `gtf2gff.pl --gff3` command for the 202 representative GFF3 files |
| `10_release_packaging/derive_subphylum_map.py` | reconstructed | subphylum map and Record No. -> label table from Workbook 3 |
| `05_representative_isoforms/build_prefilter_view.py` | new_validation | pre-filter view for scripts that read the post-filter folders |
| `10_release_packaging/check_release_inputs.py` | new_validation | read-only ID consistency check of the per-genome files |
| `tests/` | new_validation | light-weight tests; the syntax audit of `historical/` files runs only with `--include-historical` |
| `docs/USAGE_NOTES.md` | new | usage notes |
