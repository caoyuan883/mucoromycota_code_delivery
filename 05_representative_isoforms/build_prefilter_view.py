#!/usr/bin/env python3
"""Assemble a read-only *prefilter view* of the 202 representative files in a
temporary test folder: 192 unchanged final files + 10 prefilter originals.

Status: new_validation. Nothing is copied or modified: the view consists of
symbolic links (or, with --copy, copies) inside a NEW folder. The script never
writes into --final-dir or --replace-dir.

File names are matched exactly as  <label><suffix>; no prefix/fuzzy matching.
The 10 replaced labels are read from --replace-list (one label per line).

Exit codes: 0 ok; 1 output folder exists; 2 missing/extra/duplicate files.

Example (historical folder layout):
  python3 build_prefilter_view.py \
    --labels inputs/release_labels_202.txt \
    --final-dir  $REP_ROOT/02_202_longest_isoform_proteins --final-suffix _longest_isoform.faa \
    --replace-dir $REP_ROOT/05_10_prefilter/10_prefilter_faa --replace-suffix _longest_isoform.faa \
    --replace-list inputs/filtered10_labels.txt \
    --out /tmp/prefilter_view_faa
"""
import argparse, hashlib, os, shutil, sys


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_list(p):
    with open(p, encoding="utf-8") as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", required=True, help="202 release labels, one per line")
    ap.add_argument("--final-dir", required=True)
    ap.add_argument("--final-suffix", required=True, help="e.g. _longest_isoform.faa or _longest_isoform.gtf")
    ap.add_argument("--replace-dir", required=True, help="folder with the 10 prefilter originals")
    ap.add_argument("--replace-suffix", required=True)
    ap.add_argument("--replace-list", required=True, help="the 10 post-filtered labels")
    ap.add_argument("--out", required=True, help="NEW folder for the view (must not exist)")
    ap.add_argument("--out-suffix", default=None, help="name of links in the view (default: --final-suffix)")
    ap.add_argument("--copy", action="store_true", help="copy instead of symlink")
    ap.add_argument("--no-hash", action="store_true", help="skip SHA256 in the manifest")
    a = ap.parse_args()

    labels, repl = read_list(a.labels), read_list(a.replace_list)
    errs = []
    for name, lst in (("labels", labels), ("replace-list", repl)):
        d = sorted({x for x in lst if lst.count(x) > 1})
        if d:
            errs.append(f"duplicated entries in {name}: {d}")
    extra = sorted(set(repl) - set(labels))
    if extra:
        errs.append(f"replace-list labels not in labels: {extra}")
    for d in (a.final_dir, a.replace_dir):
        if not os.path.isdir(d):
            errs.append(f"folder not found: {d}")
    if errs:
        for e in errs:
            print("ERROR: " + e, file=sys.stderr)
        sys.exit(2)
    if os.path.exists(a.out):
        print(f"ERROR: output exists, refusing to overwrite: {a.out}", file=sys.stderr)
        sys.exit(1)

    rs = set(repl)
    plan = []
    for lab in labels:
        src = os.path.join(a.replace_dir, lab + a.replace_suffix) if lab in rs else os.path.join(a.final_dir, lab + a.final_suffix)
        if not os.path.isfile(src):
            errs.append(f"missing: {src}")
        plan.append((lab, "prefilter_original" if lab in rs else "final_unchanged", src))
    # files in replace-dir that are not expected
    exp_r = {lab + a.replace_suffix for lab in repl}
    unexp = sorted(f for f in os.listdir(a.replace_dir) if f.endswith(a.replace_suffix) and f not in exp_r)
    if unexp:
        errs.append(f"unexpected files in replace-dir: {unexp}")
    if errs:
        for e in errs:
            print("ERROR: " + e, file=sys.stderr)
        sys.exit(2)

    os.makedirs(a.out)
    suf = a.out_suffix or a.final_suffix
    with open(os.path.join(a.out, "00_prefilter_view_manifest.tsv"), "w", encoding="utf-8") as fh:
        fh.write("label\torigin\tsource_path\tsha256\n")
        for lab, origin, src in plan:
            dst = os.path.join(a.out, lab + suf)
            if a.copy:
                shutil.copy2(src, dst)
            else:
                os.symlink(os.path.abspath(src), dst)
            fh.write(f"{lab}\t{origin}\t{os.path.abspath(src)}\t{'-' if a.no_hash else sha256(src)}\n")
    n_r = sum(1 for p in plan if p[1] == "prefilter_original")
    print(f"view: {a.out}  files={len(plan)}  final_unchanged={len(plan) - n_r}  prefilter_original={n_r}")


if __name__ == "__main__":
    main()
