#!/usr/bin/env python3
"""Derive the strain-label -> subphylum table (00_subphylum_map.tsv) needed by
pack_release_202.sh.

Status: reconstructed. The historical 00_subphylum_map.tsv was not preserved.
This script derives it from Workbook 3 (Record No., names, Subphylum) and the
202 release labels (folder/file prefixes). It does not guess: every record must
match exactly one label and every label exactly one record, otherwise it stops.

Matching rules, tried for both "Scientific name (Workbook 1)" and
"Name used in this study":
  R1 exact      : spaces -> underscores
  R2 sp_joined  : "sp. X" -> "sp.X", then spaces -> underscores
  R3 alnum      : case-insensitive comparison after removing every
                  non-alphanumeric character
A record is ambiguous if the rules (over both name fields) produce more than one
distinct label. Optional --overrides (Record No.<TAB>label, author-confirmed)
take precedence and are reported.

Exit codes: 0 ok; 2 input error; 3 unmatched/ambiguous/non-bijective;
4 subphylum counts differ from --expect.
"""
import argparse, csv, os, re, sys
from collections import Counter, defaultdict

ABBR = {"Mucoromycotina": "Muc", "Mortierellomycotina": "Mor", "Glomeromycotina": "Glo"}


def norm_alnum(s):
    return re.sub(r"[^0-9a-z]", "", s.lower())


def candidates(name, labels, alnum_index):
    out = []
    if not name:
        return out
    r1 = name.strip().replace(" ", "_")
    if r1 in labels:
        out.append((r1, "R1_exact"))
    r2 = re.sub(r"\bsp\.\s+", "sp.", name.strip()).replace(" ", "_")
    if r2 in labels:
        out.append((r2, "R2_sp_joined"))
    for lab in alnum_index.get(norm_alnum(name), []):
        out.append((lab, "R3_alnum"))
    return out


def read_labels(args):
    if args.labels:
        with open(args.labels, encoding="utf-8") as fh:
            labs = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    else:
        if not os.path.isdir(args.labels_from_dir):
            sys.exit(f"ERROR: folder not found: {args.labels_from_dir}")
        labs = []
        for f in sorted(os.listdir(args.labels_from_dir)):
            p = os.path.join(args.labels_from_dir, f)
            if args.suffix:
                if f.endswith(args.suffix):
                    labs.append(f[: -len(args.suffix)])
            elif os.path.isdir(p) and f not in ABBR:
                labs.append(f)
    dup = [k for k, v in Counter(labs).items() if v > 1]
    if dup:
        print(f"ERROR: duplicated labels: {dup}", file=sys.stderr)
        sys.exit(2)
    return labs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workbook3", required=True, help="Workbook3_Final_202.tsv")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--labels", help="text file, one release label per line")
    g.add_argument("--labels-from-dir", help="release folder (one sub-folder per strain) or file folder with --suffix")
    ap.add_argument("--suffix", default="", help="file suffix when --labels-from-dir points to files, e.g. _longest_isoform.faa")
    ap.add_argument("--overrides", help="optional author-confirmed Record No.<TAB>label table")
    ap.add_argument("--expect", default="Mucoromycotina=123,Mortierellomycotina=57,Glomeromycotina=22")
    ap.add_argument("--out", required=True, help="output map (strain<TAB>subphylum, Muc/Mor/Glo)")
    ap.add_argument("--report", help="per-record report TSV (default: OUT with .report.tsv)")
    ap.add_argument("--record-map", help="optional output: Record No.<TAB>label<TAB>rule (record-to-label table)")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    if os.path.exists(a.out) and not a.overwrite:
        sys.exit(f"ERROR: {a.out} exists (use --overwrite)")
    labels = read_labels(a)
    lset = set(labels)
    alnum_index = defaultdict(list)
    for l in labels:
        alnum_index[norm_alnum(l)].append(l)

    with open(a.workbook3, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    need = ["Record No.", "Scientific name (Workbook 1)", "Name used in this study", "Subphylum"]
    miss = [c for c in need if rows and c not in rows[0]]
    if not rows or miss:
        sys.exit(f"ERROR: Workbook 3 missing columns {miss}")
    overrides = {}
    if a.overrides:
        with open(a.overrides, encoding="utf-8") as fh:
            for line in fh:
                if line.strip() and not line.startswith(("#", "Record")):
                    k, v = line.rstrip("\r\n").split("\t")[:2]
                    overrides[k.strip()] = v.strip()

    report, problems = [], []
    rec2lab = {}
    for r in rows:
        rec = r["Record No."].strip()
        sub = r["Subphylum"].strip()
        if sub not in ABBR:
            problems.append(f"record {rec}: unknown subphylum '{sub}'")
        if rec in overrides:
            lab = overrides[rec]
            if lab not in lset:
                problems.append(f"record {rec}: override label '{lab}' not among labels")
            rec2lab[rec] = lab
            report.append((rec, r["Scientific name (Workbook 1)"], r["Name used in this study"], sub, lab, "override", "OK"))
            continue
        c = candidates(r["Scientific name (Workbook 1)"], lset, alnum_index) + \
            candidates(r["Name used in this study"], lset, alnum_index)
        labs = sorted({x for x, _ in c})
        rules = ",".join(sorted({y for _, y in c}))
        if len(labs) == 1:
            rec2lab[rec] = labs[0]
            st = "OK"
        elif not labs:
            st = "UNMATCHED"
            problems.append(f"record {rec}: no label matched")
        else:
            st = "AMBIGUOUS"
            problems.append(f"record {rec}: ambiguous labels {labs}")
        report.append((rec, r["Scientific name (Workbook 1)"], r["Name used in this study"], sub,
                       ";".join(labs) if labs else "-", rules or "-", st))

    lab_count = Counter(rec2lab.values())
    for lab, n in lab_count.items():
        if n > 1:
            problems.append(f"label {lab} claimed by {n} records: {[k for k, v in rec2lab.items() if v == lab]}")
    unused = sorted(lset - set(rec2lab.values()))
    for lab in unused:
        problems.append(f"label {lab} matched by no record")

    rep = a.report or re.sub(r"\.tsv$", "", a.out) + ".report.tsv"
    with open(rep, "w", encoding="utf-8") as fh:
        fh.write("record_no\tscientific_name_wb1\tname_used\tsubphylum\tlabel\trule\tstatus\n")
        for t in report:
            fh.write("\t".join(t) + "\n")

    if problems:
        for p in problems:
            print("ERROR: " + p, file=sys.stderr)
        print(f"STOP: {len(problems)} problem(s); map not written. See {rep}", file=sys.stderr)
        sys.exit(3)

    sub_of = {r["Record No."].strip(): r["Subphylum"].strip() for r in rows}
    counts = Counter(sub_of[k] for k in rec2lab)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("strain\tsubphylum\n")
        for rec in sorted(rec2lab, key=lambda x: int(x) if x.isdigit() else x):
            fh.write(f"{rec2lab[rec]}\t{ABBR[sub_of[rec]]}\n")
    if a.record_map:
        with open(a.record_map, "w", encoding="utf-8") as fh:
            fh.write("Record No.\tlabel\trule\n")
            for t in report:
                fh.write(f"{t[0]}\t{t[4]}\t{t[5]}\n")
    rule_counts = Counter(t[5] for t in report)
    print(f"records={len(rows)} labels={len(labels)} mapped={len(rec2lab)}")
    print("rules: " + "; ".join(f"{k}={v}" for k, v in sorted(rule_counts.items())))
    print("counts: " + ", ".join(f"{k}={counts.get(k, 0)}" for k in ABBR))
    exp = dict(kv.split("=") for kv in a.expect.split(",")) if a.expect else {}
    bad = {k: (counts.get(k, 0), int(v)) for k, v in exp.items() if counts.get(k, 0) != int(v) or len(rec2lab) != sum(map(int, exp.values()))}
    if bad:
        print(f"ERROR: subphylum counts differ from expectation: {bad}", file=sys.stderr)
        sys.exit(4)
    print(f"written: {a.out}  report: {rep}")


if __name__ == "__main__":
    main()
