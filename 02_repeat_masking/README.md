# 02 Repeat masking

**Purpose.** Build species-specific repeat libraries (RepeatModeler2),
soft-mask the genomes (RepeatMasker `-xsmall`) and summarise repeat content.

**Inputs.** Genome FASTA files (user-supplied folder).

**Dependencies.** RepeatModeler 2.0.7, RepeatMasker 4.2.3, RMBlast 2.14.1+,
GNU parallel, Python 3 (see `docs/DEPENDENCIES.md`).

**Commands.**
```bash
bash run_repeat_softmasking.sh GENOME_DIR OUT_DIR THREADS_PER_JOB   # PARALLEL_JOBS=10 inside the script
bash summarize_repeatmodeler_logs.sh SRC_DIR OUT_CSV THREADS_PER_JOB
python3 check_softmask.py --help
```

**Order.** `run_repeat_softmasking.sh` -> `summarize_repeatmodeler_logs.sh`;
`check_softmask.py` checks the soft-masked FASTA files.

**Outputs.** `01_softmasked_genomes/`, `02_repeat_working_dir/`, run-time summary CSV.

## Repeat tables

Two table builders are kept in `historical/` with their original paths:

| entry point | inputs | outputs |
|---|---|---|
| `historical/build_table_S5_repeat_masking.py` (`--tbl-dir --rmod-dir --families-dir --log-dir --outdir`; optional `--out-dir`) | per-genome RepeatMasker `.tbl` files, RepeatModeler logs, RepeatModeler `families.fa`, pipeline logs; optional RepeatMasker `.out` files | `S5_master_table.tsv`, `S5_parse_report.txt`; with `.out` files also `S5_repeat_landscape_long.tsv` |
| `historical/S5_build/s5_v2.py` -> `fix_s5.py` -> `add_cov.py` -> `build_final_S5.py` | the same `.tbl`, RepeatModeler log, `families.fa` and pipeline-log folders; a soft-mask count table (`--softmask`) for the last step | `S5_master_table.tsv` + `S5_parse_report.txt` -> `S5_final_table.tsv` -> `S5_final_table_v2.tsv` -> `TableS5_final.tsv` |

`S5_build/cmp_mask.py A.fa B.fa` compares the soft-masked positions of two FASTA files.

The two builders write different column sets; `S5_build/` also restricts the
rows to the released genomes and adds soft-mask consistency columns.

## Script names in run logs

The run logs and `historical/` copies use the earlier names below (full list:
`docs/RENAME_MAP.tsv`). Output folders, model names and data labels keep their
original names.

| release path | historical path | note |
|---|---|---|
| `run_repeat_softmasking.sh` | `06_run_repeat_softmasking.sh` |  |
| `summarize_repeatmodeler_logs.sh` | `07_extract_rmod_log_summary.sh` |  |
