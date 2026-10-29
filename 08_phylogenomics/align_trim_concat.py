#!/usr/bin/env python3
"""
Step 02 v2 (run LOCALLY): align -> trim -> concatenate, for the P side and the
G side with byte-identical settings.

Changes vs v1
-------------
1. NEXUS charset names are prefixed with "L" (locus IDs start with a digit).
2. Resume: loci whose trimmed file already exists and is non-empty are skipped
   (use --force to recompute everything).
3. Per-locus failures are captured and logged instead of killing the whole run.
4. New report: per-taxon non-gap fraction in each final supermatrix.
5. Tool paths (not just versions) are recorded in the run log.

Environment
-----------
source /home/py/bioinfo/miniforge3/etc/profile.d/conda.sh
conda activate /home/py/bioinfo/miniforge3/envs/gene_tree      # MAFFT 7.526 + ClipKIT 2.13.2

Usage
-----
python3 align_trim_concat.py \
    --root /home/py/Desktop/Dual_tree_2026 \
    --locus-list /home/py/Desktop/Dual_tree_2026/reports/locus_list_occ80.txt \
    --tag occ80 --threads 15 --trimmer clipkit

Outputs (under <root>/03_<tag>)
-------------------------------
  align_P/ align_G/                 MAFFT output per locus
  trim_P/  trim_G/                  trimmed alignments
  supermatrix_P_<tag>.fasta         concatenated matrix, '-' for missing taxa
  supermatrix_G_<tag>.fasta
  partitions_P_<tag>.nex            kept for reviewers; NOT used in the main run
  partitions_G_<tag>.nex
  locus_stats_<tag>.tsv             per locus: n_taxa, raw len, trimmed len
  taxon_coverage_<tag>.tsv          per taxon: non-gap chars and fraction, P and G
  failed_loci_<tag>.tsv             loci that errored, with the error text
  run_log_<tag>.txt                 tool paths, versions, parameters, totals
"""

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

MAFFT_ARGS = ["--auto", "--anysymbol", "--quiet", "--thread", "1"]
CLIPKIT_ARGS = ["-m", "smart-gap"]
TRIMAL_ARGS = ["-automated1"]
MIN_TRIM_LEN = 50


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


def process_locus(job):
    locus, src, adir, tdir, trimmer, force = job
    a_out = os.path.join(adir, locus + ".aln")
    t_out = os.path.join(tdir, locus + ".trim")
    try:
        done = (not force and os.path.exists(t_out) and os.path.getsize(t_out) > 0
                and os.path.exists(a_out) and os.path.getsize(a_out) > 0)
        if not done:
            with open(a_out, "w") as fh:
                subprocess.run(["mafft"] + MAFFT_ARGS + [src], stdout=fh,
                               stderr=subprocess.DEVNULL, check=True)
            if trimmer == "clipkit":
                subprocess.run(["clipkit", a_out] + CLIPKIT_ARGS + ["-o", t_out],
                               check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.PIPE)
            else:
                subprocess.run(["trimal", "-in", a_out, "-out", t_out] + TRIMAL_ARGS,
                               check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.PIPE)
        raw = read_fasta(a_out)
        trm = read_fasta(t_out)
        raw_len = len(next(iter(raw.values()))) if raw else 0
        trim_len = len(next(iter(trm.values()))) if trm else 0
        return {"locus": locus, "ok": True, "n_taxa": len(trm),
                "raw_len": raw_len, "trim_len": trim_len, "err": ""}
    except subprocess.CalledProcessError as e:
        msg = (e.stderr.decode(errors="replace").strip().replace("\t", " ").replace("\n", " ")
               if e.stderr else "returncode=%s" % e.returncode)
        return {"locus": locus, "ok": False, "n_taxa": 0, "raw_len": 0,
                "trim_len": 0, "err": msg[:300]}
    except Exception as e:  # noqa: BLE001
        return {"locus": locus, "ok": False, "n_taxa": 0, "raw_len": 0,
                "trim_len": 0, "err": ("%s: %s" % (type(e).__name__, e))[:300]}


def run_side(side, root, out, loci, threads, trimmer, tag, force):
    src_dir = os.path.join(root, "02_paired_%s" % side)
    adir = os.path.join(out, "align_%s" % side)
    tdir = os.path.join(out, "trim_%s" % side)
    os.makedirs(adir, exist_ok=True)
    os.makedirs(tdir, exist_ok=True)

    jobs = []
    missing = []
    for locus in loci:
        src = os.path.join(src_dir, locus + ".faa")
        if os.path.exists(src):
            jobs.append((locus, src, adir, tdir, trimmer, force))
        else:
            missing.append(locus)

    stats = {}
    failed = []
    with ProcessPoolExecutor(max_workers=threads) as ex:
        for r in ex.map(process_locus, jobs):
            if r["ok"]:
                stats[r["locus"]] = (r["n_taxa"], r["raw_len"], r["trim_len"])
            else:
                failed.append((r["locus"], r["err"]))

    keep = [l for l in loci if stats.get(l, (0, 0, 0))[2] >= MIN_TRIM_LEN]
    short = [l for l in loci if l in stats and stats[l][2] < MIN_TRIM_LEN]

    taxa = set()
    for locus in keep:
        taxa |= set(read_fasta(os.path.join(tdir, locus + ".trim")))
    taxa = sorted(taxa)

    parts = []
    cat = {t: [] for t in taxa}
    pos = 1
    for locus in keep:
        seqs = read_fasta(os.path.join(tdir, locus + ".trim"))
        L = len(next(iter(seqs.values())))
        for t in taxa:
            cat[t].append(seqs.get(t, "-" * L))
        parts.append((locus, pos, pos + L - 1))
        pos += L
    matrix = {t: "".join(cat[t]) for t in taxa}
    total = pos - 1

    write_fasta(os.path.join(out, "supermatrix_%s_%s.fasta" % (side, tag)), matrix)
    with open(os.path.join(out, "partitions_%s_%s.nex" % (side, tag)), "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for locus, a, b in parts:
            fh.write("  charset L%s = %d-%d;\n" % (locus.replace(".", "_"), a, b))
        fh.write("end;\n")

    cov = {}
    gaps = 0
    for t in taxa:
        s = matrix[t]
        g = s.count("-")
        gaps += g
        cov[t] = (total - g, (total - g) / total if total else 0.0)

    return {"side": side, "loci_in": len(jobs), "loci_kept": len(keep),
            "loci_short": short, "loci_missing": missing, "failed": failed,
            "taxa": len(taxa), "taxa_list": taxa, "length": total,
            "gap_frac": gaps / (len(taxa) * total) if taxa and total else 0.0,
            "stats": stats, "keep": keep, "cov": cov}


def tool_info(cmd):
    path = shutil.which(cmd) or "not found"
    ver = "unknown"
    for args in ([cmd, "--version"], [cmd, "-v"], [cmd, "--v"]):
        try:
            r = subprocess.run(args, capture_output=True, text=True, timeout=30)
            lines = [x for x in (r.stdout + r.stderr).splitlines() if x.strip()]
            if lines:
                ver = lines[0].strip()
                break
        except Exception:  # noqa: BLE001
            continue
    return path, ver


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--locus-list", required=True)
    ap.add_argument("--tag", default="occ80")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--trimmer", choices=["clipkit", "trimal"], default="clipkit")
    ap.add_argument("--force", action="store_true", help="recompute existing alignments")
    args = ap.parse_args()

    if not shutil.which("mafft"):
        sys.exit("ERROR: mafft not found on PATH (activate the gene_tree env)")
    if not shutil.which(args.trimmer):
        sys.exit("ERROR: %s not found on PATH" % args.trimmer)

    out = os.path.join(args.root, "03_%s" % args.tag)
    os.makedirs(out, exist_ok=True)
    with open(args.locus_list) as fh:
        loci = [l.strip() for l in fh if l.strip()]

    res = {}
    for side in ("P", "G"):
        res[side] = run_side(side, args.root, out, loci, args.threads,
                             args.trimmer, args.tag, args.force)

    with open(os.path.join(out, "locus_stats_%s.tsv" % args.tag), "w") as fh:
        fh.write("locus\tn_taxa_P\traw_len_P\ttrim_len_P\tn_taxa_G\traw_len_G\ttrim_len_G\n")
        for locus in loci:
            p = res["P"]["stats"].get(locus, (0, 0, 0))
            g = res["G"]["stats"].get(locus, (0, 0, 0))
            fh.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (locus, p[0], p[1], p[2], g[0], g[1], g[2]))

    all_taxa = sorted(set(res["P"]["taxa_list"]) | set(res["G"]["taxa_list"]))
    with open(os.path.join(out, "taxon_coverage_%s.tsv" % args.tag), "w") as fh:
        fh.write("taxon\tnongap_P\tfrac_P\tnongap_G\tfrac_G\n")
        for t in all_taxa:
            p = res["P"]["cov"].get(t, (0, 0.0))
            g = res["G"]["cov"].get(t, (0, 0.0))
            fh.write("%s\t%d\t%.4f\t%d\t%.4f\n" % (t, p[0], p[1], g[0], g[1]))

    with open(os.path.join(out, "failed_loci_%s.tsv" % args.tag), "w") as fh:
        fh.write("side\tlocus\terror\n")
        for side in ("P", "G"):
            for locus, err in res[side]["failed"]:
                fh.write("%s\t%s\t%s\n" % (side, locus, err))

    mafft_p, mafft_v = tool_info("mafft")
    trim_p, trim_v = tool_info(args.trimmer)
    lines = [
        "locus list      : %s (%d loci)" % (args.locus_list, len(loci)),
        "mafft           : %s | %s | args: %s" % (mafft_p, mafft_v, " ".join(MAFFT_ARGS)),
        "trimmer         : %s | %s | args: %s" % (
            trim_p, trim_v, " ".join(CLIPKIT_ARGS if args.trimmer == "clipkit" else TRIMAL_ARGS)),
        "min trimmed len : %d" % MIN_TRIM_LEN,
        "",
    ]
    for side in ("P", "G"):
        r = res[side]
        lines.append("%s side: loci_in=%d loci_kept=%d taxa=%d length=%d gap_frac=%.4f "
                     "failed=%d too_short=%d missing_input=%d"
                     % (side, r["loci_in"], r["loci_kept"], r["taxa"], r["length"],
                        r["gap_frac"], len(r["failed"]), len(r["loci_short"]),
                        len(r["loci_missing"])))
    same_loci = res["P"]["keep"] == res["G"]["keep"]
    same_taxa = res["P"]["taxa_list"] == res["G"]["taxa_list"]
    lines += ["", "identical locus set on both sides: %s" % same_loci,
              "identical taxon set on both sides: %s" % same_taxa]
    if not same_loci or not same_taxa:
        lines.append("WARNING: sides are not symmetric -> fix before tree inference")

    low = [(t, res["P"]["cov"].get(t, (0, 0))[1], res["G"]["cov"].get(t, (0, 0))[1])
           for t in all_taxa]
    low = [x for x in low if min(x[1], x[2]) < 0.05]
    if low:
        lines.append("")
        lines.append("WARNING: taxa with <5% non-gap characters (IQ-TREE may reject):")
        for t, fp, fg in low:
            lines.append("  %-45s P=%.4f G=%.4f" % (t, fp, fg))

    txt = "\n".join(lines) + "\n"
    with open(os.path.join(out, "run_log_%s.txt" % args.tag), "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
