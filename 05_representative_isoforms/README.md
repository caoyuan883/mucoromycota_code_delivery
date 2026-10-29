# 05 Representative isoforms

**Purpose.** Keep one representative transcript per gene, write the matching
GTF and FAA, and convert the representative GTF to structural GFF3.

**Selection rule** (`filter_longest_isoforms_from_gtf.py`): total CDS length
(descending) -> transcript span (descending) -> first appearance in the GTF ->
transcript ID. `repair_3_longest_isoform_samples.py` was used for three genomes
(Linnemannia_amoeboidea, Rhizopus_microsporus, Rhizopus_stolonifer) and selects
the longest CDS with an available protein.

**Inputs.** BRAKER GTF `<Strain>.gtf` and proteins `<Strain>_protein.faa`
(from `04_structural_annotation/extract_braker_gtf_from_archives.py` or the
packaged BRAKER results).

**Dependencies.** Python 3 standard library; `gtf_to_gff3.sh` needs Perl and
AUGUSTUS `gtf2gff.pl`.

**Commands.**
```bash
python3 filter_longest_isoforms_from_gtf.py --gtf-dir GTF --protein-dir PROT --out-root OUT
python3 repair_3_longest_isoform_samples.py --gtf-dir GTF --protein-dir PROT --out-root EXISTING_OUT
G2G=.../gtf2gff.pl GTF_DIR=.../01_202_longest_isoform_gtf GFF3_DIR=NEW_DIR bash gtf_to_gff3.sh
python3 build_prefilter_view.py --help
```

`--overwrite` lets `filter_longest_isoforms_from_gtf.py` write into an existing
output folder.

**Order.** filter -> repair (3 genomes) -> stage 06 post-filtering ->
`gtf_to_gff3.sh` on the 202 final GTF files.

**GTF to GFF3.** For each GTF, `gtf_to_gff3.sh` runs
```bash
perl "$G2G" --gff3 --out="$output_gff3" < "$input_gtf"
```
without `--printExon` or `--printIntron`, and writes `00_gtf_to_gff3_log.tsv`
and `00_tool_version.txt` into the GFF3 folder. The historical converter
revision was not recorded, so output from another `gtf2gff.pl` revision may
differ in detail from the released GFF3 files.

**Pre-filter view.** Several stage 06/07 scripts read `01_202_*`/`02_202_*`,
which hold the post-filter files. `build_prefilter_view.py` builds a temporary
folder of symbolic links: the 192 unchanged final files plus the 10 pre-filter
originals (`05_10_prefilter/`), matched by exact file name and the list
`inputs/filtered10_labels.txt`. Source files are not modified.

**Format checks.** `check_gtf.py GTF_DIR FAA_DIR` and `check_faa.py FAA_DIR`
write `01_gtf_manifest.txt`/`01_gtf_check_report.txt` and
`02_manifest.tsv`/`02_check_report.txt` into the parent folder of the first
argument and overwrite files of the same name. Run them on a copy or on a view
from `build_prefilter_view.py`. For a read-only check of a release folder use
`10_release_packaging/check_release_inputs.py`.

**Outputs.** `01_longest_isoform_gtf/`, `02_longest_isoform_proteins/`,
`03_mapping_tables/`, `04_reports/`; GFF3 folder with log and version file.

**Tests.** Tie-breaks F01/F02; GFF3 wrapper G01-G04.
