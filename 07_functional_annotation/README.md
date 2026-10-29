# 07 Functional annotation

**Purpose.** Freeze the InterProScan and eggNOG-mapper results for the
representative proteins, compute per-genome statistics before and after the
post-filtering of 10 genomes, and assemble the per-genome release tables.

**Dependencies.** InterProScan 5.78-109.0 (InterPro release 109) with 8 member
analyses (`-appl Pfam,PANTHER,SUPERFAMILY,CDD,SMART,ProSitePatterns,ProSiteProfiles,AntiFam -goterms -dp -f TSV`);
eggNOG-mapper 2.1.13 (`--itype proteins -m diamond --cpu 16`, from the
annotation file headers); Python 3, bash, gzip, md5sum or md5. The InterProScan
version and member database versions are recorded in
`evidence/interproscan_VERSIONS.txt`, a text block written by
`freeze_interproscan_results.sh`, and in the installation path used by the run
scripts. See `docs/DEPENDENCIES.md`.

## Order

1. Freeze and repair (`ips_fix_20261011/`; env `W`, `IPS`; `FREEZE_ROOT` or `F`):
   `freeze_emapper_results.sh` -> `diagnose_interproscan.sh` -> `rerun_interproscan_split.sh`
   (copies and calls `run_chunk.sh` from its own folder; keep the folder layout) ->
   `check_interproscan_results.sh` -> `freeze_interproscan_results.sh`.
   Without `F`, the frozen folder is named `<FREEZE_ROOT>/{emapper,interproscan}_v1_<today>`;
   set `F` to reuse a fixed name such as `interproscan_v1_20261011`.
2. Statistics:
   `summarize_functional_annotation.py --emapper-dir E --ips-dir I --faa PREFILTER_FAA --out OUT`
   (pre-filter, 202 genomes; exits 1 unless every genome has both results), then
   `recompute_filtered10_functional_stats.py --emapper-dir E --ips-dir I --filtered FILTERED_FAA --table OUT/01_per_strain_functional_summary.tsv --orig PREFILTER_FAA --out OUT2`
   -> `02_per_strain_functional_summary_final.tsv` (copy in `inputs/`).
3. Release tables: `FUNC_ROOT=... FAADIR=FINAL_FAA bash assemble_functional_tables.sh`;
   `FUNC_ROOT=... FILTDIR=FILTERED_FAA bash subset_filtered10_functional_tables.sh`
   (`STRICT=1`: exit 3 on any problem).
4. Independent recheck:
   `recheck_func_numbers.py --ips-dir IPS --em-dir EM --prefilter-faa-dir PREFILTER_FAA --filtered-dir FILTERED_FAA --out OUT`
   (optional: `--subphylum-map`, `--fail-on-diff`).
   `--prefilter-faa-dir` must hold the pre-filter proteins; `02_202_*` holds
   the post-filter set, so use a view from
   `05_representative_isoforms/build_prefilter_view.py`. Without
   `--subphylum-map` the subphylum is assigned from a genus list; with a map
   that does not cover all genomes the script exits 4.

**Definitions.** Denominator = representative proteins of the release FAA;
hit = at least one line in the tool output; union and intersection over protein
IDs; GO terms by source (`InterPro`, `PANTHER`) for InterProScan and the `GOs`
column for eggNOG-mapper. Tests A01 (synthetic) and A02 (set identities on the
final table).

**Subphylum medians.** `recheck_func_numbers.py` writes `subphylum_medians.tsv`
with the per-subphylum median of the union coverage (`pct_any`, percent of
proteins with an InterProScan or eggNOG-mapper hit) over the released genomes.
The median is `statistics.median` (for an even number of genomes, the mean of
the two middle values); the table gives the unrounded value and the value
rounded to one decimal (half up). For Glomeromycotina (22 genomes) the
unrounded median is 72.252888 (mean of 71.761577 and 72.744199), shown as
72.3%; Mucoromycotina 91.2%, Mortierellomycotina 82.4%. The `EXPECT` block of
the script now holds 72.3 for `median_any_Glo` (previously 72.2); the value is
still recomputed and compared at one decimal. `evidence/recheck_report.txt` is
the historical report written with the earlier value and therefore lists
`median_any_Glo` as DIFF (72.2 vs 72.25). `summarize_functional_annotation.py` reports
cohort-wide statistics with a nearest-rank median in `05_cohort_summary.tsv`.

**Outputs.** Frozen tool outputs, `01_`/`02_per_strain_functional_summary*.tsv`,
per-genome release tables, `subphylum_medians.tsv`, `evidence/recheck_report.txt`,
`evidence/07_build_log.tsv`.

## Script names in run logs

The run logs and `historical/` copies use the earlier names below (full list:
`docs/RENAME_MAP.tsv`). Output folders, model names and data labels keep their
original names.

| release path | historical path | note |
|---|---|---|
| `assemble_functional_tables.sh` | `make_07_func.sh` |  |
| `ips_fix_20261011/check_interproscan_results.sh` | `ips_fix_20261011/fix_20261011_part2/04_check_ips.sh` |  |
| `ips_fix_20261011/diagnose_interproscan.sh` | `ips_fix_20261011/02_diagnose_ips.sh` |  |
| `ips_fix_20261011/freeze_emapper_results.sh` | `ips_fix_20261011/01_freeze_emapper.sh` |  |
| `ips_fix_20261011/freeze_interproscan_results.sh` | `ips_fix_20261011/fix_20261011_part2/06_freeze_ips.sh` |  |
| `ips_fix_20261011/rerun_interproscan_split.sh` | `ips_fix_20261011/03_rerun_ips_split.sh` |  |
| `recompute_filtered10_functional_stats.py` | `func_filtered10.py` |  |
| `subset_filtered10_functional_tables.sh` | `make_filtered10_func.sh` |  |
| `summarize_functional_annotation.py` | `func_summary_v3.py` |  |
