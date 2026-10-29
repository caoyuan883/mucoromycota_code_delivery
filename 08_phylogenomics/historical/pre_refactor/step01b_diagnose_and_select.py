#!/usr/bin/env python3
"""
Step 01b: diagnose the QC log, profile every taxon, and write locus lists
for the occupancy thresholds you may use in step 02.

Reads only what step01 already produced. Nothing is modified or deleted.

Usage
-----
python3 step01b_diagnose_and_select.py --out /home/py/Desktop/Dual_tree_2026

Outputs (under <out>/reports)
-----------------------------
  qc_reason_counts.tsv     drop reasons x side (multi_record rows separated out)
  taxon_stats.tsv          per-taxon: loci in P, in G, paired, occupancy, drop counts
  locus_list_occ50.txt / occ70 / occ80 / occ90 / occ95   locus IDs per threshold
  threshold_report.tsv     for each threshold: #loci, #taxa below 50% coverage, est. aa
  step01b_summary.txt      the printed report
"""

import argparse
import os
import statistics
from collections import Counter, defaultdict

THRESHOLDS = [0.50, 0.70, 0.80, 0.90, 0.95]


def read_matrix(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        taxa = header[1:]
        mat = {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            mat[parts[0]] = {taxa[i]: parts[1 + i] == "1" for i in range(len(taxa))}
    return mat, taxa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rep = os.path.join(args.out, "reports")

    P, taxa = read_matrix(os.path.join(rep, "presence_P.tsv"))
    G, _ = read_matrix(os.path.join(rep, "presence_G.tsv"))
    loci = sorted(set(P) & set(G))
    n_tax = len(taxa)

    paired = {l: {t for t in taxa if P[l].get(t) and G[l].get(t)} for l in loci}

    reason_counts = Counter()
    taxon_drops = defaultdict(Counter)
    with open(os.path.join(rep, "qc_dropped.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            side, taxon, locus, reason = f[0], f[1], f[2], f[3]
            key = reason.split("_median")[0] if reason.startswith("length_outlier") else reason
            reason_counts[(side, key)] += 1
            if key != "multi_record_longest_kept":
                taxon_drops[taxon][side] += 1

    with open(os.path.join(rep, "qc_reason_counts.tsv"), "w") as fh:
        fh.write("side\treason\tcount\n")
        for (side, reason), n in sorted(reason_counts.items(), key=lambda x: -x[1]):
            fh.write("%s\t%s\t%d\n" % (side, reason, n))

    occ_by_locus = {l: len(paired[l]) / n_tax for l in loci}
    med_len = {}
    occ_path = os.path.join(rep, "occupancy_paired.tsv")
    if os.path.exists(occ_path):
        with open(occ_path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5:
                    med_len[f[0]] = int(f[3])

    rows = []
    for t in taxa:
        nP = sum(1 for l in loci if P[l].get(t))
        nG = sum(1 for l in loci if G[l].get(t))
        nPair = sum(1 for l in loci if t in paired[l])
        rows.append((t, nP, nG, nPair, nPair / len(loci), taxon_drops[t]["P"], taxon_drops[t]["G"]))
    rows.sort(key=lambda r: r[4])
    with open(os.path.join(rep, "taxon_stats.tsv"), "w") as fh:
        fh.write("taxon\tloci_P\tloci_G\tloci_paired\tfrac_paired\tdropped_P\tdropped_G\n")
        for r in rows:
            fh.write("%s\t%d\t%d\t%d\t%.4f\t%d\t%d\n" % r)

    lines = ["=== drop reasons (side, reason, count) ==="]
    for (side, reason), n in sorted(reason_counts.items(), key=lambda x: -x[1]):
        lines.append("%s\t%s\t%d" % (side, reason, n))

    lines.append("")
    lines.append("=== 15 taxa with lowest paired-locus coverage ===")
    for r in rows[:15]:
        lines.append("%-45s paired=%4d  frac=%.3f  dropP=%4d dropG=%4d" % (r[0], r[3], r[4], r[5], r[6]))

    lines.append("")
    lines.append("=== threshold report ===")
    trows = []
    for th in THRESHOLDS:
        keep = [l for l in loci if occ_by_locus[l] >= th]
        with open(os.path.join(rep, "locus_list_occ%d.txt" % int(th * 100)), "w") as fh:
            fh.write("\n".join(keep) + "\n")
        if keep:
            per_tax = [sum(1 for l in keep if t in paired[l]) / len(keep) for t in taxa]
            low = sum(1 for x in per_tax if x < 0.50)
            est = sum(med_len.get(l, 0) for l in keep)
            trows.append((th, len(keep), min(per_tax), statistics.median(per_tax), low, est))
            lines.append(
                "occ>=%.2f  loci=%4d  taxon_coverage min=%.3f median=%.3f  taxa<50%%=%d  est_aa~%d"
                % (th, len(keep), min(per_tax), statistics.median(per_tax), low, est)
            )
    with open(os.path.join(rep, "threshold_report.tsv"), "w") as fh:
        fh.write("threshold\tn_loci\tmin_taxon_cov\tmedian_taxon_cov\ttaxa_below_50pct\test_aa_before_trim\n")
        for r in trows:
            fh.write("%.2f\t%d\t%.4f\t%.4f\t%d\t%d\n" % r)

    txt = "\n".join(lines) + "\n"
    with open(os.path.join(rep, "step01b_summary.txt"), "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
