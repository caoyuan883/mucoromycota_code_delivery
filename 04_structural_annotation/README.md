# 04 Structural annotation (BRAKER2 and BRAKER3)

**Purpose.** Gene prediction on soft-masked genomes with two independent batch
systems, kept in separate folders with their own configuration:

* `braker2/` - BRAKER2 with protein evidence (genomes without usable RNA-seq);
* `braker3/` - BRAKER3 with RNA-seq (merged BAM from stage 03) and protein evidence.

The system used for each genome is assigned in stage 03 (`map_species_braker_mode.sh`).

**Inputs.**
* Soft-masked genomes below `DATA_ROOT`: `03_softmasked_genomes_v3/` (BRAKER2) and
  `03_softmasked_genomes/` (BRAKER3).
* Protein evidence `Fungi.clean.fa`, prepared once from OrthoDB v12 Fungi `Fungi.fa`
  by cutting each header to its first field:
  `awk '/^>/{print $1; next} {print}' Fungi.fa > Fungi.clean.fa`.
* Merged RNA-seq BAM per species (BRAKER3) and the species list
  (`species_mode_mapping.csv` or `SPECIES_MODE_CSV`).
* Local `mucoromycota_odb12` lineage for compleasm and best_by_compleasm.

Genome FASTA headers are cleaned per genome by `clean_genome_headers.sh`.

**Dependencies.** BRAKER, AUGUSTUS, GeneMark (ES/EP/ETP), ProtHint, DIAMOND,
compleasm, TSEBRA with a modified `best_by_compleasm.py` (see below), SAMtools.
Versions: `docs/DEPENDENCIES.md`.

**Configuration.** Both `config.sh` files require `ENV_ROOT`, `DATA_ROOT` and
`BUSCO_LINEAGE_DIR`; export them or set `BRAKER_SITE_ENV` to an edited copy of
`config/braker_site.env.example`. Concurrency is set by `THREADS_PER_JOB` and
`MAX_PARALLEL` (example setting: 8 threads per job, up to 4 jobs in parallel, as
in the usual runs; use `MAX_PARALLEL=1` when resources are limited or for single
restarts). Further variables: `BBC_CWD` (`06b_bbc_merge.sh`, `06c_bbc_full.sh`;
folder holding the compleasm `mb_downloads`, default: script folder),
`BATCH_ROOT` (`monitor_addon_two.sh`, `monitor_addon_one.sh`) and
`SPECIES_MODE_CSV` (`braker3/prepare_species_list.sh`, default
`$DATA_ROOT/species_mode_mapping.csv`).

## BRAKER2 (`braker2/`)

`config.sh` (sourced by the other scripts) -> `clean_genome_headers.sh` ->
`preflight.sh` -> `run_batch.sh` (runs `run_braker2_one.sh` per genome) ->
`finalize_manifest.sh` (writes `_manifest_final.tsv`) -> `package_results.sh` ->
`make_raw_gff3.sh` -> `write_package_docs.sh`.

Add-on runs: `addon_two_genomes.sh` (monitor `monitor_addon_two.sh`) and
`addon_one_genome.sh` (monitor `monitor_addon_one.sh`). Housekeeping:
`archive_intermediates.sh`, `stage_cleanup.sh`, `archive_addon_two.sh`.

best_by_compleasm helpers (historical names): `06_bbc.sh` and `06_bbc_all.sh`
score existing protein sets with compleasm; `06b_bbc_merge.sh` and
`06c_bbc_full.sh` back up `braker.gtf`, `braker.aa` and `bbc/`
(`*.backup_pre_merge`) and run `best_by_compleasm.py` with output in
`_bbc_merge/` or `_bbc_full/`, without copying results back to
`braker.gtf`/`braker.aa`; `09_rebuild_manifest.sh` writes `_manifest_v4.tsv`.

## BRAKER3 (`braker3/`)

`config.sh` (sourced) -> `prepare_species_list.sh` -> `clean_genome_headers.sh` ->
`preflight.sh` -> `clean_augustus_config.sh` -> `run_batch.sh` (runs
`run_braker3_one.sh`; a run that ends as `DONE_PARTIAL` exits 6).
`evidence/braker3_bbc_summary_2026-09-25.tsv` summarises the best_by_compleasm
log lines per genome, including the re-runs made with
`historical/10_rerun_bbc_serial.sh`.

## best_by_compleasm post-processing

BRAKER gene sets were post-processed with `best_by_compleasm.py` (TSEBRA) and the
`mucoromycota_odb12` lineage. The script was modified locally to skip HMM hits
without a score cutoff, and a `lengths_cutoff` file missing from the lineage was
generated locally (HMM ID, 0.8 x HMM length, 1, 0). Because the HMM names of this
lineage did not match the identifiers of its score-cutoff table, no BUSCO
transcripts were passed on: the step that adds BUSCO transcripts to the BRAKER
gene set added none. The re-merging step, which enforces the AUGUSTUS or GeneMark
gene set when the BRAKER set lacks more than 20% of BUSCOs and is at least five
points worse (default settings), still ran.

The re-merged gene sets were adopted for *Cunninghamella echinulata*,
*Mucor heterogamus*, *Rhizopus microsporus*, *Rhizopus stolonifer* and
*Linnemannia amoeboidea*. For *Umbelopsis ramanniana* a re-merged set was
produced, but the BRAKER2 re-annotation was used.

Running this step requires the locally modified `best_by_compleasm.py` (production
revision: skip HMM hits without a cutoff, SHA256
`b51574065bc255e424a2eadf713b7905dbbe2e1c168177293e9c99f929c99f4d`). It is not
included in this repository or in standard BRAKER installations; see
`docs/DEPENDENCIES.md`.

## After both systems

`extract_braker_gtf_from_archives.py --source-dirs BRAKER2_CORE BRAKER3_DELIVERY --out OUT`
collects `braker.gtf` from the packaged results for stage 05.

**Outputs.** Per genome `braker.aa`, `braker.gtf`, `braker.gff3`.

`braker2/make_raw_gff3.sh` runs `gtf2gff.pl --gff3 --printExon --printIntron` to
write BRAKER raw GFF3 files. The GFF3 files of the 202 representative gene sets
are written by `05_representative_isoforms/gtf_to_gff3.sh`.

**Historical copies.** `historical/` contains earlier batch and packaging
scripts, the superseded packager `braker2_batch_v3/07_package_core.sh` and
`braker3_env.sh` (not used by the batch scripts).

## Script names in run logs

The run logs and `historical/` copies use the earlier names below (full list:
`docs/RENAME_MAP.tsv`). Output folders, model names and data labels keep their
original names.

| release path | historical path | note |
|---|---|---|
| `braker2/06_bbc.sh` | `braker2_batch_v3/06_bbc.sh` | historical name kept |
| `braker2/06_bbc_all.sh` | `braker2_batch_v3/06_bbc_all.sh` | historical name kept |
| `braker2/06b_bbc_merge.sh` | `braker2_batch_v3/06b_bbc_merge.sh` | historical name kept |
| `braker2/06c_bbc_full.sh` | `braker2_batch_v3/06c_bbc_full.sh` | historical name kept |
| `braker2/09_rebuild_manifest.sh` | `braker2_batch_v3/09_rebuild_manifest.sh` | historical name kept |
| `braker2/addon_one_genome.sh` | `braker2_batch_v3/17_addon_one.sh` |  |
| `braker2/addon_two_genomes.sh` | `braker2_batch_v3/14_addon_two.sh` |  |
| `braker2/archive_addon_two.sh` | `braker2_batch_v3/16_archive_two_new.sh` |  |
| `braker2/archive_intermediates.sh` | `braker2_batch_v3/12_archive_selective.sh` |  |
| `braker2/clean_genome_headers.sh` | `braker2_batch_v3/01_clean_headers.sh` |  |
| `braker2/config.sh` | `braker2_batch_v3/00_config.sh` |  |
| `braker2/finalize_manifest.sh` | `braker2_batch_v3/09b_finalize_manifest.sh` |  |
| `braker2/make_raw_gff3.sh` | `braker2_batch_v3/10_make_gff3.sh` |  |
| `braker2/monitor_addon_one.sh` | `braker2_batch_v3/18_watch_one.sh` |  |
| `braker2/monitor_addon_two.sh` | `braker2_batch_v3/15_watch_addon.sh` |  |
| `braker2/package_results.sh` | `braker2_batch_v3/07_package_core_v4.sh` |  |
| `braker2/preflight.sh` | `braker2_batch_v3/02_preflight.sh` |  |
| `braker2/run_batch.sh` | `braker2_batch_v3/04_supervisor.sh` |  |
| `braker2/run_braker2_one.sh` | `braker2_batch_v3/03_run_braker2_one.sh` |  |
| `braker2/stage_cleanup.sh` | `braker2_batch_v3/13_stage_cleanup.sh` |  |
| `braker2/write_package_docs.sh` | `braker2_batch_v3/11_toplevel_docs.sh` |  |
| `braker3/clean_augustus_config.sh` | `braker3_batch_v4/02_clean_augustus_config.sh` |  |
| `braker3/clean_genome_headers.sh` | `braker3_batch_v4/00b_clean_headers.sh` |  |
| `braker3/config.sh` | `braker3_batch_v4/00_config.sh` |  |
| `braker3/evidence/braker3_bbc_summary_2026-09-25.tsv` | `braker3_batch_v4/evidence/braker3_bbc_summary_2026-09-25.tsv` |  |
| `braker3/preflight.sh` | `braker3_batch_v4/01_preflight.sh` |  |
| `braker3/prepare_species_list.sh` | `braker3_batch_v4/00_prep_species_list.sh` |  |
| `braker3/run_batch.sh` | `braker3_batch_v4/04_supervisor.sh` |  |
| `braker3/run_braker3_one.sh` | `braker3_batch_v4/03_run_braker3_one.sh` |  |
| `historical/braker2_batch_v3/07_package_core.sh` | `braker2_batch_v3/07_package_core.sh` | moved to historical, unchanged |
