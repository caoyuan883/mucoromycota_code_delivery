#!/usr/bin/env python3
"""
Step 01 for the dual-tree (protein-mode vs genome-mode) BUSCO validation.

What it does
------------
1. Scans two BUSCO single-copy directories laid out as  <root>/<Strain>/<BUSCOID>.faa
2. Renames every sequence header to the strain (folder) name  -> Newick-safe leaf label
3. Runs identical sequence QC on both sides (internal stop codons, X fraction,
   length outliers, multi-record files)
4. Writes per-locus merged FASTA for each side (all sequences that passed QC)
5. Applies the PAIRWISE-VALID rule: a (locus, taxon) unit is kept only if it
   passed QC on BOTH sides; otherwise it is recorded as missing on both sides
6. Writes occupancy tables so you can choose the occupancy threshold in step 02

No third-party dependencies (pure standard library).

Usage
-----
python3 collect_busco_markers.py \
    --p-dir /home/py/Desktop/Mucoromycota_protein_fungi_odb12_single_copy_BUSCO \
    --g-dir /home/py/Desktop/Mucoromycota_genome_fungi_odb12_single_copy_BUSCO \
    --out   /home/py/Desktop/Dual_tree_2026

Outputs (under --out)
---------------------
  01_merged_P/<BUSCOID>.faa        all QC-passing P sequences, headers = leaf labels
  01_merged_G/<BUSCOID>.faa        same for G
  02_paired_P/<BUSCOID>.faa        pairwise-valid subset, common loci only  <-- use these
  02_paired_G/<BUSCOID>.faa        same for G
  reports/taxa_P.txt, taxa_G.txt, taxa_mismatch.txt
  reports/qc_dropped.tsv           every dropped sequence + reason
  reports/presence_P.tsv           locus x taxon 0/1 matrix (post-QC)
  reports/presence_G.tsv
  reports/occupancy_paired.tsv     per-locus taxon count in the paired set
  reports/summary.txt              numbers to quote in the manuscript
"""

import argparse
import os
import statistics
import sys
from collections import defaultdict

NEWICK_BAD = " \t()[]{}':;,|"


def safe_label(name):
    out = []
    for ch in name:
        out.append("_" if ch in NEWICK_BAD else ch)
    return "".join(out)


def read_fasta(path):
    recs = []
    header = None
    chunks = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    recs.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        recs.append((header, "".join(chunks)))
    return recs


def qc_sequence(seq, max_x_frac):
    """Return (clean_seq, reason_or_None). reason != None means: drop."""
    s = seq.strip().upper().replace("-", "").replace(".", "")
    while s.endswith("*"):
        s = s[:-1]
    if not s:
        return None, "empty"
    if "*" in s:
        return None, "internal_stop"
    n_x = s.count("X") + s.count("*") + s.count("?")
    if len(s) < 30:
        return None, "too_short_lt30aa"
    if n_x / len(s) > max_x_frac:
        return None, "x_fraction_gt_%.2f" % max_x_frac
    return s, None


def scan_side(root, tag, max_x_frac, dropped_rows):
    """-> data[locus][taxon] = seq ; taxa set"""
    data = defaultdict(dict)
    taxa = set()
    if not os.path.isdir(root):
        sys.exit("ERROR: directory not found: %s" % root)
    for entry in sorted(os.listdir(root)):
        sdir = os.path.join(root, entry)
        if not os.path.isdir(sdir):
            continue
        label = safe_label(entry)
        taxa.add(label)
        for fn in sorted(os.listdir(sdir)):
            if not fn.endswith((".faa", ".fa", ".fasta")):
                continue
            locus = fn.rsplit(".", 1)[0]
            recs = read_fasta(os.path.join(sdir, fn))
            if not recs:
                dropped_rows.append([tag, label, locus, "no_record", ""])
                continue
            if len(recs) > 1:
                recs = sorted(recs, key=lambda r: len(r[1]), reverse=True)
                dropped_rows.append(
                    [tag, label, locus, "multi_record_longest_kept", str(len(recs))]
                )
            orig_id, raw = recs[0]
            clean, reason = qc_sequence(raw, max_x_frac)
            if reason:
                dropped_rows.append([tag, label, locus, reason, orig_id])
                continue
            if label in data[locus]:
                dropped_rows.append([tag, label, locus, "duplicate_label", orig_id])
                continue
            data[locus][label] = clean
    return data, taxa


def flag_length_outliers(data, tag, dropped_rows, lo=0.5, hi=2.0):
    """Drop sequences whose length deviates strongly from the locus median."""
    for locus in list(data.keys()):
        lens = [len(s) for s in data[locus].values()]
        if len(lens) < 4:
            continue
        med = statistics.median(lens)
        if med <= 0:
            continue
        for taxon in list(data[locus].keys()):
            L = len(data[locus][taxon])
            if L < lo * med or L > hi * med:
                dropped_rows.append(
                    [tag, taxon, locus, "length_outlier_%d_vs_median_%d" % (L, med), ""]
                )
                del data[locus][taxon]


def write_matrix(path, data, taxa):
    loci = sorted(data.keys())
    taxa = sorted(taxa)
    with open(path, "w") as fh:
        fh.write("locus\t" + "\t".join(taxa) + "\n")
        for locus in loci:
            row = ["1" if t in data[locus] else "0" for t in taxa]
            fh.write(locus + "\t" + "\t".join(row) + "\n")


def dump_fasta_dir(outdir, data, loci=None, wrap=60):
    os.makedirs(outdir, exist_ok=True)
    for locus in sorted(loci if loci is not None else data.keys()):
        seqs = data.get(locus, {})
        if not seqs:
            continue
        with open(os.path.join(outdir, locus + ".faa"), "w") as fh:
            for taxon in sorted(seqs):
                s = seqs[taxon]
                fh.write(">" + taxon + "\n")
                for i in range(0, len(s), wrap):
                    fh.write(s[i:i + wrap] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p-dir", required=True)
    ap.add_argument("--g-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-x-frac", type=float, default=0.05)
    args = ap.parse_args()

    rep = os.path.join(args.out, "reports")
    os.makedirs(rep, exist_ok=True)
    dropped = []

    P, taxaP = scan_side(args.p_dir, "P", args.max_x_frac, dropped)
    G, taxaG = scan_side(args.g_dir, "G", args.max_x_frac, dropped)
    flag_length_outliers(P, "P", dropped)
    flag_length_outliers(G, "G", dropped)

    with open(os.path.join(rep, "taxa_P.txt"), "w") as fh:
        fh.write("\n".join(sorted(taxaP)) + "\n")
    with open(os.path.join(rep, "taxa_G.txt"), "w") as fh:
        fh.write("\n".join(sorted(taxaG)) + "\n")
    only_p = sorted(taxaP - taxaG)
    only_g = sorted(taxaG - taxaP)
    with open(os.path.join(rep, "taxa_mismatch.txt"), "w") as fh:
        fh.write("# leaf labels present in only one side -> MUST be fixed before step 02\n")
        for t in only_p:
            fh.write("P_only\t%s\n" % t)
        for t in only_g:
            fh.write("G_only\t%s\n" % t)

    write_matrix(os.path.join(rep, "presence_P.tsv"), P, taxaP)
    write_matrix(os.path.join(rep, "presence_G.tsv"), G, taxaG)
    dump_fasta_dir(os.path.join(args.out, "01_merged_P"), P)
    dump_fasta_dir(os.path.join(args.out, "01_merged_G"), G)

    common_taxa = taxaP & taxaG
    common_loci = sorted(set(P.keys()) & set(G.keys()))

    pairedP = defaultdict(dict)
    pairedG = defaultdict(dict)
    for locus in common_loci:
        shared = (set(P[locus]) & set(G[locus])) & common_taxa
        for taxon in shared:
            pairedP[locus][taxon] = P[locus][taxon]
            pairedG[locus][taxon] = G[locus][taxon]

    kept_loci = [l for l in common_loci if pairedP[l]]
    dump_fasta_dir(os.path.join(args.out, "02_paired_P"), pairedP, kept_loci)
    dump_fasta_dir(os.path.join(args.out, "02_paired_G"), pairedG, kept_loci)

    n_tax = len(common_taxa)
    with open(os.path.join(rep, "occupancy_paired.tsv"), "w") as fh:
        fh.write("locus\tn_taxa\toccupancy\tmedian_len_P\tmedian_len_G\n")
        for locus in kept_loci:
            n = len(pairedP[locus])
            mp = statistics.median([len(s) for s in pairedP[locus].values()])
            mg = statistics.median([len(s) for s in pairedG[locus].values()])
            fh.write("%s\t%d\t%.4f\t%d\t%d\n" % (locus, n, n / n_tax if n_tax else 0, mp, mg))

    with open(os.path.join(rep, "qc_dropped.tsv"), "w") as fh:
        fh.write("side\ttaxon\tlocus\treason\tnote\n")
        for row in dropped:
            fh.write("\t".join(row) + "\n")

    occ = [len(pairedP[l]) / n_tax for l in kept_loci] if n_tax else []
    lines = [
        "taxa in P: %d" % len(taxaP),
        "taxa in G: %d" % len(taxaG),
        "taxa shared: %d" % n_tax,
        "taxa only in P: %s" % (", ".join(only_p) or "none"),
        "taxa only in G: %s" % (", ".join(only_g) or "none"),
        "loci in P: %d" % len(P),
        "loci in G: %d" % len(G),
        "loci in P and G: %d" % len(common_loci),
        "loci retained after pairwise rule: %d" % len(kept_loci),
        "dropped sequence records (all reasons): %d" % len(dropped),
    ]
    if occ:
        lines += [
            "paired occupancy min/median/max: %.3f / %.3f / %.3f"
            % (min(occ), statistics.median(occ), max(occ)),
            "loci with occupancy == 1.00 : %d" % sum(1 for x in occ if x >= 0.9999),
            "loci with occupancy >= 0.90 : %d" % sum(1 for x in occ if x >= 0.90),
            "loci with occupancy >= 0.80 : %d" % sum(1 for x in occ if x >= 0.80),
            "loci with occupancy >= 0.70 : %d" % sum(1 for x in occ if x >= 0.70),
            "loci with occupancy >= 0.50 : %d" % sum(1 for x in occ if x >= 0.50),
        ]
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(rep, "summary.txt"), "w") as fh:
        fh.write(txt)
    sys.stdout.write(txt)


if __name__ == "__main__":
    main()
