# 10 Release assembly and packaging

**Purpose.** Assemble the five files per genome, check them, and pack one
`tar.gz` per genome into the three subphylum folders of the released dataset.

**Dependencies.** bash, coreutils (`md5sum` or macOS `md5`), tar, gzip, Python 3.

## Order

1. Subphylum map:
   ```bash
   python3 10_release_packaging/derive_subphylum_map.py \
       --workbook3 inputs/Workbook3_Final_202.tsv --labels inputs/release_labels_202.txt \
       --out 00_subphylum_map.tsv --record-map record_label_map.tsv
   ```
   Matching rules: exact name, `sp. X` -> `sp.X`, alphanumeric-only comparison,
   applied to both name columns of Workbook 3; one label per record and one
   record per label, otherwise exit 3; counts must be 123/57/22, otherwise
   exit 4. With the current Workbook 3 export all 202 records match one label
   each. Records 41, 94, 111, 127, 151 and 152 use the source name (naming
   exception, `docs/DATA_LINKAGE.md`); their labels are also listed in
   `inputs/label_overrides_6.tsv`. Adding `--overrides inputs/label_overrides_6.tsv`
   gives the same map (test M06).
2. Read-only check of the per-genome files:
   `python3 10_release_packaging/check_release_inputs.py --labels inputs/release_labels_202.txt --faa-dir ... --gtf-dir ... --gff3-dir ... --expect-total 2618164 --out report.tsv`
   (protein IDs FAA<->GTF and FAA<->GFF3 in both directions, duplicate IDs,
   empty sequences, CRLF, .gz; exit 5 on failure).
3. Assemble: `REP_ROOT=... FUNC_RELEASE=... OUT=NEW_RELEASE STRICT=1 bash 10_release_packaging/make_release_202.sh`
   writes the five files per genome, `00_release_check.tsv` and
   `MD5SUMS_files.txt`; exit 1 if `OUT` exists; with `STRICT=1`, exit 3 unless
   202 genomes / 2,618,164 proteins / 0 mismatches (without `STRICT=1` the
   totals are printed but not enforced).
4. Pack: `RELEASE=NEW_RELEASE SUBPHYLUM_MAP=00_subphylum_map.tsv STRICT=1 bash 10_release_packaging/pack_release_202.sh`
   writes `Mucoromycotina/`, `Mortierellomycotina/` and `Glomeromycotina/` with
   one `<label>.tar.gz` per genome, `MD5SUMS.txt` and `00_pack_log.tsv`.
   Exit 2 if the release folder or map is missing, 1 if archives exist and
   `OVERWRITE=1` is not set; genomes without a map entry are logged as
   `NO_SUBPHYLUM`; with `STRICT=1`, exit 3 unless 202 archives are OK with 123/57/22.

## Checksum lists

| file | written by | content |
|---|---|---|
| `MD5SUMS_files.txt` | `make_release_202.sh` | MD5 of each of the five files per genome |
| `MD5SUMS.txt` | `pack_release_202.sh` | MD5 of each per-genome `tar.gz` archive |

Re-packing writes new timestamps into the archives, so archive checksums of a
new run differ from those of the released archives; compare the per-file
checksums instead.

**Evidence.** `evidence/software_versions.txt` (tool versions collected in the
original environment).
