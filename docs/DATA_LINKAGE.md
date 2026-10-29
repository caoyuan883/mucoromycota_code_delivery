# Data linkage

## Identifiers

| identifier | where | example | rule |
|---|---|---|---|
| Record No. | Workbook 3 and later workbooks | `8` | 1-202, the record key of the dataset |
| Source identifier | Workbook 3 | `GCA_037043795.1`, JGI portal ID | links to the public assembly; never changed |
| release label | file and folder names | `Acaulospora_colombiana` | one per record; used in every per-genome file name |
| subphylum | Workbook 3 `Subphylum` | Mucoromycotina | 123 / 57 / 22 |
| gene / transcript ID | inside GTF, GFF3, FAA, TSV | `g1`, `g1.t1` | unique only within one genome; combine with the label |

Workbook 3 defines which records are released and the names used in this
study. Some historical scripts retain batch-specific sample counts. The released
dataset comprises the 202 genome records listed in Workbook 3.

## Record -> label -> subphylum

`10_release_packaging/derive_subphylum_map.py` links each Workbook 3 record to
exactly one release label (rules and exit codes in
`10_release_packaging/README.md`) and writes the label -> subphylum table used
for packaging; `--record-map` also writes the record -> label table. The mapping
is never inferred from row order, protein counts or name similarity.

Naming exception: six records (41, 94, 111, 127, 151, 152) keep the source
name as the name used in this study and as their release label (Workbook 3,
sheet `Naming_notes_6`, exported as `inputs/Workbook3_Naming_notes_6.tsv`). In
the current Workbook 3 export these records match their labels through the
normal rules; `inputs/label_overrides_6.tsv` lists the same six labels and can
be passed with `--overrides` as an explicit cross-check.

## Per-genome files

```
<Subphylum>/<label>.tar.gz
  <label>/<label>_longest_isoform.gtf             representative transcripts (structural)
  <label>/<label>_longest_isoform.gff3            converted from the GTF, structural only
  <label>/<label>_longest_isoform.faa             representative proteins (one per gene)
  <label>/<label>_interproscan.tsv                InterProScan 5.78-109.0 TSV (8 member analyses)
  <label>/<label>_emapper.annotations.tsv         eggNOG-mapper 2.1.13 annotations
```

Join keys within one genome:

* Structural files: the FAA header ID, the GTF `transcript_id` and the GFF3 mRNA
  `ID` denote the same representative transcript and its protein. Both
  `10_release_packaging/check_release_inputs.py` and `make_release_202.sh`
  compare these ID sets.
* Functional tables contain only proteins with a reported result. In the
  InterProScan TSV, column 1 is the protein ID and one protein can have several
  rows (one per match); in the eggNOG-mapper table, `query` is the protein ID.
  These ID sets are subsets of the FAA IDs, not equal to them.
  `07_functional_annotation/assemble_functional_tables.sh` counts functional
  IDs that are absent from the FAA (expected 0); the release check scripts
  above do not read the functional tables.

## Pre-filter and release

* Release: 202 genomes, 2,618,164 proteins (192 unchanged + 10 post-filtered).
* Pre-filter: 2,639,091 proteins (the same 202 genomes before 20,927 proteins
  were removed from 10 genomes; labels in `inputs/filtered10_labels.txt`).
* Statistics that refer to the pre-filter state use a pre-filter view
  (`05_representative_isoforms/build_prefilter_view.py`).
* The phylogenomic matrices were built from the BUSCO runs described in
  `08_phylogenomics/README.md`.

## Workbooks and scripts

| workbook | content | produced or checked by |
|---|---|---|
| Workbook 3 | genome records, source identifiers, names, subphylum | curated table; exports in `inputs/Workbook3_Final_202.tsv` and `inputs/Workbook3_Naming_notes_6.tsv` (source in `inputs/README.md`) |
| Workbook 5 | repeat content | `02_repeat_masking/historical/` table builders (see `02_repeat_masking/README.md`) |
| Workbook 7 | pre-filter annotation structure and the changes in the 10 filtered genomes | `06_quality_and_filtering/workbook7_sheet2_summarize_filtered_annotations.py` |
| Workbook 8 | final BUSCO / OMArk quality | `06_quality_and_filtering/` |
| Workbook 9 | functional annotation per genome | `07_functional_annotation/recompute_filtered10_functional_stats.py` (output copy: `inputs/02_per_strain_functional_summary_final.tsv`) |

## Small inputs in `inputs/`

See `inputs/README.md`.
