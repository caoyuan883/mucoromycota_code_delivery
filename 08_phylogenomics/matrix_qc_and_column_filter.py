#!/usr/bin/env python3
"""
Step 02c: diagnose the gap composition of each supermatrix and (optionally)
apply a symmetric column filter.

Why: a raw "gap fraction" mixes two very different things.
  * structural gaps  = the taxon has no sequence for that locus at all
                       (an entire locus block is gaps)  -> real missing data
  * alignment gaps   = the taxon HAS a sequence, but the aligner inserted gaps
                       -> a property of the alignment, not of the data
Only the first is "missing data" in the sense of the phylogenomics literature.

The column filter removes alignment columns whose gap fraction exceeds
--max-col-gap (default 0.90). The SAME rule is applied to the P and G matrices.

Usage
-----
python3 matrix_qc_and_column_filter.py --root /home/py/Desktop/Dual_tree_2026 \
        --tag occ80 --max-col-gap 0.9 --extra-thresholds 0.7,0.5

Outputs (under <root>/03_<tag>)
-------------------------------
  matrix_qc_<tag>.txt                     the printed report
  column_gap_histogram_<tag>.tsv          gap-fraction distribution of columns
  supermatrix_<side>_<tag>_g<NN>.fasta    filtered matrices (one per threshold)
  partitions_<side>_<tag>_g<NN>.nex       remapped partition files
"""

import argparse
import os


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs, wrap=60):
    with open(path, "w") as fh:
        for k in sorted(seqs):
            fh.write(">" + k + "\n")
            s = seqs[k]
            for i in range(0, len(s), wrap):
                fh.write(s[i:i + wrap] + "\n")


def read_partitions(path):
    blocks = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("charset"):
                name, rng = line[len("charset"):].split("=")
                a, b = rng.strip().rstrip(";").split("-")
                blocks.append((name.strip(), int(a), int(b)))
    return blocks


def analyse(seqs, blocks):
    taxa = sorted(seqs)
    n = len(taxa)
    L = len(seqs[taxa[0]])
    col_gaps = [0] * L
    for t in taxa:
        s = seqs[t]
        for i, ch in enumerate(s):
            if ch == "-":
                col_gaps[i] += 1
    total_gaps = sum(col_gaps)

    structural = 0
    absent_units = 0
    total_units = 0
    for _, a, b in blocks:
        span = b - a + 1
        for t in taxa:
            total_units += 1
            if seqs[t].count("-", a - 1, b) == span:
                structural += span
                absent_units += 1
    align_gaps = total_gaps - structural
    return {
        "taxa": n, "length": L, "col_gaps": col_gaps,
        "gap_frac": total_gaps / (n * L),
        "structural_frac": structural / (n * L),
        "align_frac": align_gaps / (n * L),
        "unit_absent_frac": absent_units / total_units if total_units else 0.0,
    }


def filter_columns(seqs, blocks, max_col_gap):
    taxa = sorted(seqs)
    n = len(taxa)
    L = len(seqs[taxa[0]])
    col_gaps = [0] * L
    for t in taxa:
        s = seqs[t]
        for i, ch in enumerate(s):
            if ch == "-":
                col_gaps[i] += 1
    keep = [i for i in range(L) if col_gaps[i] / n <= max_col_gap]
    keep_set = set(keep)
    out = {t: "".join(seqs[t][i] for i in keep) for t in taxa}

    newblocks = []
    pos = 1
    for name, a, b in blocks:
        kept = sum(1 for i in range(a - 1, b) if i in keep_set)
        if kept > 0:
            newblocks.append((name, pos, pos + kept - 1))
            pos += kept
    return out, newblocks


def write_partitions(path, blocks):
    with open(path, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for name, a, b in blocks:
            fh.write("  charset %s = %d-%d;\n" % (name, a, b))
        fh.write("end;\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--tag", default="occ80")
    ap.add_argument("--max-col-gap", type=float, default=0.9)
    ap.add_argument("--extra-thresholds", default="",
                    help="comma-separated extra thresholds, e.g. 0.7,0.5")
    args = ap.parse_args()

    out = os.path.join(args.root, "03_%s" % args.tag)
    thresholds = [args.max_col_gap]
    if args.extra_thresholds.strip():
        thresholds += [float(x) for x in args.extra_thresholds.split(",") if x.strip()]

    lines = []
    hist_rows = []
    for side in ("P", "G"):
        mpath = os.path.join(out, "supermatrix_%s_%s.fasta" % (side, args.tag))
        ppath = os.path.join(out, "partitions_%s_%s.nex" % (side, args.tag))
        seqs = read_fasta(mpath)
        blocks = read_partitions(ppath)
        a = analyse(seqs, blocks)
        n = a["taxa"]

        lines.append("=== %s side (%s) ===" % (side, args.tag))
        lines.append("taxa                     : %d" % n)
        lines.append("columns                  : %d" % a["length"])
        lines.append("total gap fraction       : %.4f" % a["gap_frac"])
        lines.append("  structural (no seq)    : %.4f" % a["structural_frac"])
        lines.append("  alignment-inserted     : %.4f" % a["align_frac"])
        lines.append("absent taxon x locus     : %.4f" % a["unit_absent_frac"])

        bins = [0] * 11
        for g in a["col_gaps"]:
            bins[min(10, int((g / n) * 10))] += 1
        for i, c in enumerate(bins):
            hist_rows.append((side, "%.1f-%.1f" % (i / 10, (i + 1) / 10), c,
                              c / a["length"]))

        for th in thresholds:
            filt, nb = filter_columns(seqs, blocks, th)
            L2 = len(next(iter(filt.values())))
            g2 = sum(s.count("-") for s in filt.values()) / (n * L2) if L2 else 0
            cov = [1 - (s.count("-") / L2) for s in filt.values()] if L2 else [0]
            tagn = int(round(th * 100))
            write_fasta(os.path.join(out, "supermatrix_%s_%s_g%d.fasta"
                                     % (side, args.tag, tagn)), filt)
            write_partitions(os.path.join(out, "partitions_%s_%s_g%d.nex"
                                          % (side, args.tag, tagn)), nb)
            lines.append("  filter col_gap<=%.2f : columns=%d (%.1f%% kept)  "
                         "gap_frac=%.4f  min_taxon_cov=%.4f  loci_kept=%d"
                         % (th, L2, 100.0 * L2 / a["length"], g2, min(cov), len(nb)))
        lines.append("")

    with open(os.path.join(out, "column_gap_histogram_%s.tsv" % args.tag), "w") as fh:
        fh.write("side\tgap_frac_bin\tn_columns\tprop_columns\n")
        for r in hist_rows:
            fh.write("%s\t%s\t%d\t%.4f\n" % r)

    txt = "\n".join(lines) + "\n"
    with open(os.path.join(out, "matrix_qc_%s.txt" % args.tag), "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
