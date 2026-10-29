# Tests

Minimal run (Python standard library only):

```bash
python3 tests/run_tests.py --outdir /tmp/mucoro_tests
```

Full run with external inputs and optional tools:

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install dendropy==5.1.0
python3 tests/run_tests.py --outdir /tmp/mucoro_tests \
    --trees /path/to/folder/with/the/four/treefiles \
    --gtf2gff /path/to/Augustus/scripts/gtf2gff.pl \
    --rscript /path/to/Rscript \
    --plot-python /path/to/python3/with/pandas/and/matplotlib
```

Historical syntax audit (adds S03 and S04h):

```bash
python3 tests/run_tests.py --outdir /tmp/mucoro_tests_hist --include-historical
```

`--outdir` must be a new folder. Results: `OUTDIR/test_results.tsv` and
`OUTDIR/test_results.md` (command, input, expected, actual, exit code, status).
Status is `PASS`, `FAIL` or `NOT_RUN`; `NOT_RUN` is never counted as passed.
The suite exits 1 if any test fails.

The default run checks the public (non-historical) entry points. Files under
`historical/` are kept byte-identical and are checked only with
`--include-historical`: S03 (Python and shell) and S04h (R). S03 then reports
`01_data_selection/historical/07_build_species_list.py`, a historical file with
a `.py` suffix that contains shell commands; the file is kept unchanged.

Optional inputs and the tests that need them:

| missing | tests reported as `NOT_RUN` |
|---|---|
| DendroPy | T01-T03, X15 and H11 (`compare_four_trees.py --help`; the script imports DendroPy at start) |
| `--trees` | T01-T03, X15 |
| `--gtf2gff` | G03 |
| `--rscript` (and no `Rscript` on `PATH`) | S04, the R `--help` tests and the R plotting draft tests |
| `--plot-python` (and no Python with pandas and matplotlib) | the Python plotting tests |

D01-D07 (plotting drafts that need input tables not in the repository) and
N01-N11 (real-data regressions) are always `NOT_RUN`; their commands and
expected results are listed in the result table.

Reference run (Linux; Python 3.14 with DendroPy 5.1.0; R 4.4.3; system Python 3
with pandas and matplotlib; `--trees`, `--gtf2gff`, `--rscript`,
`--plot-python`): 112 tests, PASS 93, FAIL 0, NOT_RUN 19. With
`--include-historical` in addition: 113 tests, PASS 94, FAIL 1 (S03), NOT_RUN 18.

| group | what is checked | data |
|---|---|---|
| S | `ast.parse` / `bash -n` / R parse on every public script (S01, S02, S04); historical audit S03, S04h only with `--include-historical` | repository |
| H | `--help` of every argparse entry point and of the plotting scripts | - |
| C | shell entry points stop with a non-zero exit when required variables are missing | - |
| F | longest-isoform key and tie-breaks (F01-F02); pre-filter view (F03-F05) | synthetic (`SYN_`) |
| G | GTF -> GFF3 wrapper (command, no extra options, overwrite refusal) | synthetic |
| R | release file checker (ID mismatch, duplicates, empty sequences, CRLF, .gz) | synthetic |
| M | subphylum map: rules and stop conditions incl. an ambiguous record (M01-M04, synthetic); Workbook 3 export without overrides (M05), same result with `inputs/label_overrides_6.tsv` (M06), agreement with `per_strain_release.tsv` (M07) | synthetic + `inputs/` |
| A | functional statistics: denominators, union/intersection, GO split (A01-A02); subphylum medians of union coverage from `inputs/per_strain_release.tsv` (A03) and on synthetic data (A04) | synthetic + `inputs/` |
| P | `make_release_202.sh` / `pack_release_202.sh` on a synthetic release | synthetic |
| T | four tree files: tip sets and RF distances after outgroup pruning | tree files (`--trees`) |
| X | consistency: `$here`/`$SD` calls, paths in docs and config, script names, `docs/RENAME_MAP.tsv`, scripts named in READMEs, figure index and required fields, plot SHA256, Illustrator statement, plotting drafts without input (X10-X11), missing-column stops (X12-X14), caption constant 197/204 / nRF (X15), best_by_compleasm scripts (X16), relative links in public Markdown (X17), `config/phylogenomics_logged_runs.tsv` against the IQ-TREE logs and `run_iqtree_occ80.sh` (X18) | repository, tree files |
| D | plotting drafts that need input tables | `NOT_RUN`, inputs and commands listed |
| N | real-data regressions (BRAKER, BUSCO, OMArk, InterProScan, eggNOG-mapper, IQ-TREE, full release) | `NOT_RUN`, commands listed |

`synthetic_data.py` writes all synthetic inputs; every label starts with
`SYN_`. Expected values are fixed in `run_tests.py`.
