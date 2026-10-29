#!/usr/bin/env python3
"""Light-weight test suite for the release scripts.

Runs only safe checks: syntax, --help, missing-configuration exits and
SYNTHETIC inputs (tests/synthetic_data.py), plus read-only checks on the small
real inputs in inputs/ and (optionally) the four published treefiles.
It never downloads databases and never starts BRAKER, BUSCO, OMArk,
InterProScan, eggNOG-mapper or IQ-TREE. Plot drafts are written only into OUTDIR.

Usage (minimal):
  python3 tests/run_tests.py --outdir /tmp/mucoro_tests
Usage (all options):
  python3 tests/run_tests.py --outdir /tmp/mucoro_tests --trees /path/to/folder/with/4/treefiles \
      --gtf2gff /path/to/gtf2gff.pl --rscript /path/to/Rscript --plot-python /path/to/python3 \
      --include-historical
The default run checks the public entry points; --include-historical adds the
syntax audit of the files under historical/ (S03).

Writes OUTDIR/test_results.tsv and OUTDIR/test_results.md.
Status: PASS / FAIL / NOT_RUN (NOT_RUN is never counted as PASS).
Exit code: 0 if no FAIL, 1 otherwise.
"""
import argparse, csv, glob, gzip, os, platform, shutil, subprocess, sys, tarfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
R = []


def run(cmd, env=None, cwd=None, shell=False):
    e = dict(os.environ)
    if env:
        for k, v in env.items():
            if v is None:
                e.pop(k, None)
            else:
                e[k] = v
    p = subprocess.run(cmd, cwd=cwd or REPO, env=e, shell=shell, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr)


def rec(tid, cat, cmd, inp, expected, status, actual, rc):
    if not isinstance(cmd, str):
        cmd = " ".join(cmd)
    R.append(dict(id=tid, category=cat, command=cmd.replace(REPO + "/", "").replace(REPO, "."),
                  input=inp, expected=expected, actual=str(actual).replace("\n", " | ")[:400], exit_code=rc, status=status))


def case(tid, cat, cmd, inp, expected_rc, expected_txt, check=None, env=None, shell=False, cwd=None):
    rc, out = run(cmd, env=env, shell=shell, cwd=cwd)
    ok = (rc == expected_rc) if expected_rc != "nonzero" else (rc != 0)
    detail = out.strip().splitlines()[-1] if out.strip() else ""
    if ok and check:
        ok, detail = check(out)
    rec(tid, cat, cmd, inp, f"exit {expected_rc}; {expected_txt}", "PASS" if ok else "FAIL", detail, rc)
    return rc, out


def copy_plain(src, dst):
    """copy *.treefile without permissions (inputs may be read-only)"""
    os.makedirs(dst)
    for f in glob.glob(os.path.join(src, "*.treefile")):
        shutil.copyfile(f, os.path.join(dst, os.path.basename(f)))


def release_files(kind):
    out = []
    for p in sorted(glob.glob(os.path.join(REPO, "**", "*." + kind), recursive=True)):
        rel = os.path.relpath(p, REPO)
        if rel.startswith("tests/") and kind == "py":
            out.append(rel); continue
        out.append(rel)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--trees", help="folder with P_occ80_g90.treefile ... G_occ90_g90.treefile")
    ap.add_argument("--gtf2gff", help="path to AUGUSTUS gtf2gff.pl for the GFF3 conversion test")
    ap.add_argument("--rscript", help="Rscript for R parse/--help/plot-draft checks (default: Rscript on PATH)")
    ap.add_argument("--plot-python", help="python with pandas+matplotlib for the Python plotting script (default: auto)")
    ap.add_argument("--include-historical", action="store_true",
                    help="also run the syntax audit of historical/ files (S03); not part of the default run")
    a = ap.parse_args()
    if os.path.exists(a.outdir):
        sys.exit(f"ERROR: outdir exists: {a.outdir}")
    O = os.path.abspath(a.outdir); os.makedirs(O)
    S = os.path.join(O, "synthetic"); W = os.path.join(O, "work"); os.makedirs(W)
    rc, out = run([PY, "tests/synthetic_data.py", S])
    if rc:
        sys.exit("cannot create synthetic data:\n" + out)
    I = os.path.join(REPO, "inputs")

    # ---------------- syntax --------------------------------------------------
    try:
        import dendropy  # noqa
        have_dp = True
    except ImportError:
        have_dp = False
    pub = lambda fs: [f for f in fs if "/historical/" not in f]
    hist = lambda fs: [f for f in fs if "/historical/" in f]
    def py_bad(files):
        return [f for f in files if run([PY, "-c", f"import ast,sys; ast.parse(open({f!r},encoding='utf-8').read(), {f!r})"])[0]]
    files = pub(release_files("py")); bad = py_bad(files)
    rec("S01", "syntax", "python3 ast.parse on every public (non-historical) *.py", f"{len(files)} files", "all parse",
        "PASS" if not bad else "FAIL", f"{len(files) - len(bad)}/{len(files)} parse; failed: {bad}", len(bad))
    files = pub(release_files("sh")); bad = [f for f in files if run(["bash", "-n", f])[0]]
    rec("S02", "syntax", "bash -n on every public (non-historical) *.sh", f"{len(files)} files",
        "all parse", "PASS" if not bad else "FAIL", f"{len(files) - len(bad)}/{len(files)} parse; failed: {bad}", len(bad))
    # S03 historical audit: only with --include-historical (historical files are kept byte-identical)
    hpy, hsh = hist(release_files("py")), hist(release_files("sh"))
    if a.include_historical:
        bpy = py_bad(hpy); bsh = [f for f in hsh if run(["bash", "-n", f])[0]]
        rec("S03", "historical_audit", "python3 ast.parse on historical/*.py; bash -n on historical/*.sh (R files: S04h)",
            f"{len(hpy)} py + {len(hsh)} sh files", "all parse; files whose content does not match the suffix are reported",
            "PASS" if not (bpy or bsh) else "FAIL",
            f"py {len(hpy) - len(bpy)}/{len(hpy)} parse, not Python: {bpy}; sh {len(hsh) - len(bsh)}/{len(hsh)} parse, failed: {bsh}",
            len(bpy) + len(bsh))
    else:
        rec("S03", "historical_audit", "python3 ast.parse on historical/*.py; bash -n on historical/*.sh", f"{len(hpy)} py + {len(hsh)} sh files",
            "all parse", "NOT_RUN", "historical audit not requested (--include-historical)", "-")

    # ---------------- --help --------------------------------------------------
    helps = ["05_representative_isoforms/filter_longest_isoforms_from_gtf.py", "05_representative_isoforms/repair_3_longest_isoform_samples.py",
             "05_representative_isoforms/build_prefilter_view.py", "04_structural_annotation/extract_braker_gtf_from_archives.py",
             "06_quality_and_filtering/omark_filtering/audit_omark_contigs.py", "06_quality_and_filtering/omark_filtering/filter_omark_contigs.py",
             "06_quality_and_filtering/make_protein_qc_master_table.py", "07_functional_annotation/summarize_functional_annotation.py",
             "07_functional_annotation/recompute_filtered10_functional_stats.py", "07_functional_annotation/recheck_func_numbers.py",
             "08_phylogenomics/compare_four_trees.py", None,
             "10_release_packaging/derive_subphylum_map.py", "10_release_packaging/check_release_inputs.py",
             "08_phylogenomics/diagnose_and_select_markers.py", "08_phylogenomics/align_trim_concat.py", "08_phylogenomics/matrix_qc_and_column_filter.py",
             "08_phylogenomics/collect_busco_markers.py"]
    for i, h in enumerate(helps, 1):
        if h is None:  # H12 (Workbook 9 helper) moved out of the repository; ID kept unused
            continue
        if h == "08_phylogenomics/compare_four_trees.py" and not have_dp:
            rec(f"H{i:02d}", "help", f"python3 {h} --help", "-", "exit 0", "NOT_RUN", "DendroPy not installed (imported at start of script)", "-"); continue
        if not os.path.exists(os.path.join(REPO, h)):
            rec(f"H{i:02d}", "help", f"python3 {h} --help", "-", "exit 0", "FAIL", "file missing", "-"); continue
        case(f"H{i:02d}", "help", [PY, h, "--help"], "-", 0, "usage printed",
             check=lambda o: ("usage:" in o, o.strip().splitlines()[0]))

    # ---------------- missing configuration -> non-zero ----------------------
    envless = {"ENV_ROOT": None, "DATA_ROOT": None, "BUSCO_LINEAGE_DIR": None, "BRAKER_SITE_ENV": None,
               "REP_ROOT": None, "FUNC_RELEASE": None, "OUT": None, "RELEASE": None, "FUNC_ROOT": None, "FAADIR": None,
               "FILTDIR": None, "IQ": None, "WORK": None, "G2G": None, "GTF_DIR": None, "GFF3_DIR": None,
               "W": None, "IPS": None, "BATCH_ROOT": None, "FREEZE_ROOT": None}
    for i, s in enumerate(["10_release_packaging/make_release_202.sh", "10_release_packaging/pack_release_202.sh",
                           "07_functional_annotation/assemble_functional_tables.sh", "07_functional_annotation/subset_filtered10_functional_tables.sh",
                           "05_representative_isoforms/gtf_to_gff3.sh", "08_phylogenomics/run_iqtree_occ80.sh",
                           "08_phylogenomics/compare_substitution_models.sh", "04_structural_annotation/braker2/config.sh",
                           "04_structural_annotation/braker3/config.sh", "04_structural_annotation/braker2/monitor_addon_one.sh"], 1):
        case(f"C{i:02d}", "config", ["bash", os.path.join(REPO, s)], "required variables unset", "nonzero",
             "stops with a 'set ...' message, no work done", env=envless, cwd=W,
             check=lambda o: ("set " in o or "ERROR" in o, o.strip().splitlines()[-1] if o.strip() else ""))

    # ---------------- 05 longest isoform (frozen key) -------------------------
    out_li = os.path.join(W, "li")
    case("F01", "synthetic", [PY, "05_representative_isoforms/filter_longest_isoforms_from_gtf.py", "--gtf-dir", f"{S}/longest/gtf",
         "--protein-dir", f"{S}/longest/protein", "--out-root", out_li], "SYNTHETIC SYN_A.gtf: 3 genes x 2 isoforms", 0,
         "selected g1.t2 (span tie-break), g2.t1 (first appearance), g3.t2 (CDS length)",
         check=lambda o: (lambda ids: (ids == ["g1.t2", "g2.t1", "g3.t2"], f"selected {ids}"))(
             [l[1:].strip() for l in open(f"{out_li}/02_longest_isoform_proteins/SYN_A_longest_isoform.faa") if l.startswith(">")]))
    case("F02", "synthetic", [PY, "05_representative_isoforms/filter_longest_isoforms_from_gtf.py", "--gtf-dir", f"{S}/longest/gtf",
         "--protein-dir", f"{S}/longest/protein", "--out-root", out_li], "same as F01, output exists", "nonzero", "refuses non-empty output without --overwrite")

    # ---------------- 05 prefilter view ---------------------------------------
    pv = os.path.join(W, "pv")
    base = [PY, "05_representative_isoforms/build_prefilter_view.py", "--labels", f"{S}/prefilter/labels.txt", "--final-dir", f"{S}/prefilter/final",
            "--final-suffix", "_longest_isoform.faa", "--replace-dir", f"{S}/prefilter/pre", "--replace-suffix", "_longest_isoform.faa"]
    case("F03", "synthetic", base + ["--replace-list", f"{S}/prefilter/replace.txt", "--out", pv], "SYNTHETIC 3 labels, SYN_B replaced", 0,
         "SYN_B links to prefilter original, SYN_A/SYN_C to final files; sources untouched",
         check=lambda o: (open(f"{pv}/SYN_B_longest_isoform.faa").read().endswith("MPREFILTER\n") and os.path.islink(f"{pv}/SYN_A_longest_isoform.faa")
                          and open(f"{S}/prefilter/final/SYN_B_longest_isoform.faa").read().endswith("MFINALSYN_B\n"), o.strip()))
    case("F04", "synthetic", base + ["--replace-list", f"{S}/prefilter/replace.txt", "--out", pv], "same, output exists", 1, "refuses to overwrite")
    case("F05", "synthetic", base + ["--replace-list", f"{S}/prefilter/replace_missing.txt", "--out", pv + "2"],
         "replace-list names SYN_C, no SYN_C prefilter file", 2, "missing file reported (+ unexpected SYN_B file)")

    # ---------------- 05 GTF -> GFF3 wrapper ----------------------------------
    stub = os.path.join(W, "stub_gtf2gff.pl")
    with open(stub, "w") as fh:
        fh.write('# SYNTHETIC stub: records its arguments\nmy ($o)=grep {/^--out=/} @ARGV; $o=~s/^--out=//; open(my $f, ">", $o); print $f "##gff-version 3\\n# ARGS @ARGV\\n"; while(<STDIN>){}; close $f;\n')
    gdir = os.path.join(W, "gff3_stub")
    case("G01", "synthetic", ["bash", "05_representative_isoforms/gtf_to_gff3.sh"], "SYNTHETIC GTF + argument-recording stub", 0,
         "stub called exactly with --gff3 --out=... (no --printExon/--printIntron)",
         env={"G2G": stub, "GTF_DIR": f"{S}/release_check/ok/gtf", "GFF3_DIR": gdir},
         check=lambda o: (lambda t: ("# ARGS --gff3 --out=" in t and "print" not in t.split("# ARGS")[1], t.split("\n")[1]))(
             open(f"{gdir}/SYN_A_longest_isoform.gff3").read()))
    case("G02", "synthetic", ["bash", "05_representative_isoforms/gtf_to_gff3.sh"], "same, output folder not empty", 1, "refuses",
         env={"G2G": stub, "GTF_DIR": f"{S}/release_check/ok/gtf", "GFF3_DIR": gdir})
    if a.gtf2gff and os.path.isfile(a.gtf2gff):
        g3 = os.path.join(W, "gff3_real")
        case("G03", "synthetic", ["bash", "05_representative_isoforms/gtf_to_gff3.sh"], f"SYNTHETIC GTF + {a.gtf2gff}", 0, "GFF3 written",
             env={"G2G": a.gtf2gff, "GTF_DIR": f"{S}/release_check/ok/gtf", "GFF3_DIR": g3})
        case("G04", "synthetic", [PY, "10_release_packaging/check_release_inputs.py", "--labels", f"{S}/release_check/labels.txt",
             "--faa-dir", f"{S}/release_check/ok/faa", "--gtf-dir", f"{S}/release_check/ok/gtf", "--gff3-dir", g3, "--out", f"{W}/g04.tsv"],
             "SYNTHETIC FAA/GTF + GFF3 from real gtf2gff.pl", 0, "mRNA IDs of converted GFF3 equal FAA/GTF IDs both ways")
    else:
        rec("G03", "synthetic", "bash 05_representative_isoforms/gtf_to_gff3.sh with real gtf2gff.pl", "-", "exit 0", "NOT_RUN", "no --gtf2gff given", "-")

    # ---------------- 10 release triplet checker ------------------------------
    def chk(case_dir, gz=False):
        return [PY, "10_release_packaging/check_release_inputs.py", "--labels", f"{S}/release_check/labels.txt",
                "--faa-dir", f"{S}/release_check/{case_dir}/faa", "--gtf-dir", f"{S}/release_check/{case_dir}/gtf",
                "--gff3-dir", f"{S}/release_check/{case_dir}/gff3", "--out", f"{W}/chk_{case_dir}.tsv"]
    case("R01", "synthetic", chk("ok"), "SYNTHETIC consistent triplet", 0, "status OK")
    case("R02", "synthetic", chk("ok_gz"), "SYNTHETIC consistent triplet, .gz", 0, "status OK (gz read)")
    case("R03", "synthetic", chk("crlf"), "SYNTHETIC FAA with CRLF", 5, "crlf=faa reported",
         check=lambda o: ("\tfaa\t" in open(f"{W}/chk_crlf.tsv").read(), "crlf column = faa"))
    case("R04", "synthetic", chk("dup_id"), "SYNTHETIC duplicated FASTA ID", 5, "dup_faa_ids=1",
         check=lambda o: (open(f"{W}/chk_dup_id.tsv").read().splitlines()[1].split("\t")[4] == "1", "dup_faa_ids=1"))
    case("R05", "synthetic", chk("gff3_extra"), "SYNTHETIC GFF3 with extra mRNA", 5, "gff3_not_in_faa=1",
         check=lambda o: (open(f"{W}/chk_gff3_extra.tsv").read().splitlines()[1].split("\t")[9] == "1", "gff3_not_in_faa=1"))
    case("R06", "synthetic", chk("missing"), "SYNTHETIC missing GFF3", 5, "missing=gff3",
         check=lambda o: ("\tgff3\tFAIL" in open(f"{W}/chk_missing.tsv").read(), "missing=gff3"))

    # ---------------- 10 subphylum map ----------------------------------------
    sm = [PY, "10_release_packaging/derive_subphylum_map.py"]
    case("M01", "synthetic", sm + ["--workbook3", f"{S}/subphylum/wb3_ok.tsv", "--labels", f"{S}/subphylum/labels_ok.txt",
         "--expect", "Mucoromycotina=1,Mortierellomycotina=1,Glomeromycotina=1", "--out", f"{W}/m01.tsv"],
         "SYNTHETIC 3 records (exact, 'sp.X', 'KO_9'->'_KO9')", 0, "3 mapped: Muc/Mor/Glo",
         check=lambda o: (open(f"{W}/m01.tsv").read() == "strain\tsubphylum\nSynthomyces_alpha\tMuc\nSynthomyces_sp.AB-1\tMor\nSynthomyces_sp._KO9\tGlo\n", "map content as expected"))
    case("M02", "synthetic", sm + ["--workbook3", f"{S}/subphylum/wb3_ambiguous.tsv", "--labels", f"{S}/subphylum/labels_ambiguous.txt",
         "--expect", "", "--out", f"{W}/m02.tsv"], "SYNTHETIC record with two name fields -> two labels", 3, "stops, no map written",
         check=lambda o: (not os.path.exists(f"{W}/m02.tsv") and "ambiguous" in o, "ambiguous; map absent"))
    case("M03", "synthetic", sm + ["--workbook3", f"{S}/subphylum/wb3_ok.tsv", "--labels", f"{S}/subphylum/labels_dup.txt", "--out", f"{W}/m03.tsv"],
         "SYNTHETIC duplicated label list", 2, "stops")
    case("M04", "synthetic", sm + ["--workbook3", f"{S}/subphylum/wb3_ok.tsv", "--labels", f"{S}/subphylum/labels_ok.txt", "--out", f"{W}/m04.tsv"],
         "SYNTHETIC 3 records vs default expectation 123/57/22", 4, "count mismatch -> exit 4")
    EXC = {"41": "Cunninghamella_bertholletiae", "94": "Mortierella_hygrophila", "111": "Mucor_ambiguus",
           "127": "Mucor_velutinosus", "151": "Rhizopus_azygosporus", "152": "Rhizopus_delemar"}
    def recmap(path):
        return {r["Record No."]: r["label"] for r in csv.DictReader(open(path), delimiter="\t")}
    def m05check(o):
        if not os.path.exists(f"{W}/m05_records.tsv"): return False, "record map not written"
        rm = recmap(f"{W}/m05_records.tsv"); exc = {k: rm.get(k) for k in EXC}
        ok = ("Mucoromycotina=123, Mortierellomycotina=57, Glomeromycotina=22" in o and len(rm) == 202
              and len(set(rm.values())) == 202 and exc == EXC)
        return ok, f"{len(rm)} records -> {len(set(rm.values()))} labels; naming-exception records: {exc}"
    case("M05", "real_small", sm + ["--workbook3", f"{I}/Workbook3_Final_202.tsv", "--labels", f"{I}/release_labels_202.txt",
         "--out", f"{W}/m05.tsv", "--record-map", f"{W}/m05_records.tsv"],
         "inputs/Workbook3_Final_202.tsv + 202 labels, without --overrides", 0,
         "202 records -> 202 labels, 123/57/22; records 41,94,111,127,151,152 -> source-name labels", check=m05check)
    def m06check(o):
        if not (os.path.exists(f"{W}/m06.tsv") and os.path.exists(f"{W}/m05.tsv")): return False, "map missing"
        same_map = open(f"{W}/m06.tsv").read() == open(f"{W}/m05.tsv").read()
        same_rec = recmap(f"{W}/m06_records.tsv") == recmap(f"{W}/m05_records.tsv")
        return ("Mucoromycotina=123, Mortierellomycotina=57, Glomeromycotina=22" in o and same_map and same_rec,
                f"subphylum map identical to M05: {same_map}; record->label identical to M05: {same_rec}")
    case("M06", "real_small", sm + ["--workbook3", f"{I}/Workbook3_Final_202.tsv", "--labels", f"{I}/release_labels_202.txt",
         "--overrides", f"{I}/label_overrides_6.tsv", "--out", f"{W}/m06.tsv", "--record-map", f"{W}/m06_records.tsv"],
         "real inputs + inputs/label_overrides_6.tsv (cross-check)", 0, "202 mapped, 123/57/22; same map and record labels as M05", check=m06check)

    def agree(_):
        m = {r["strain"]: r["subphylum"] for r in csv.DictReader(open(f"{W}/m06.tsv"), delimiter="\t")}
        p = {r["strain"]: r["subphylum"] for r in csv.DictReader(open(f"{I}/per_strain_release.tsv"), delimiter="\t")}
        d = [k for k in p if m.get(k) != p[k]]
        return (len(m) == 202 and not d, f"{len(m)} labels; disagreements with genus-list assignment in per_strain_release.tsv: {len(d)}")
    if os.path.exists(f"{W}/m06.tsv"):
        rc_, _ = 0, ""
        ok, det = agree(None)
        rec("M07", "real_small", "compare M06 map with per_strain_release.tsv subphylum column", "M06 output", "0 disagreements", "PASS" if ok else "FAIL", det, 0)

    # ---------------- 07 functional statistics (frozen definitions) -----------
    ff = os.path.join(W, "ff")
    def ffcheck(_):
        r = {x["strain"]: x for x in csv.DictReader(open(f"{ff}/02_per_strain_functional_summary_final.tsv"), delimiter="\t")}
        a_ = r["SYN_A"]
        exp = dict(n_proteins="4", ips_hit="2", egg_hit="2", union_hit="3", both_tools="1", ips_only="1", egg_only="1",
                   go_union="2", ips_go="2", ips_go_InterPro="1", ips_go_PANTHER="1", egg_go="1", egg_ko="1", ips_Pfam="1",
                   union_pct="75.0", filtered="yes")
        bad = {k: (a_[k], v) for k, v in exp.items() if a_[k] != v}
        ok = not bad and r["SYN_B"]["filtered"] == "no" and r["SYN_B"]["n_proteins"] == "7"
        return ok, ("all expected values" if ok else f"mismatch {bad}")
    case("A01", "synthetic", [PY, "07_functional_annotation/recompute_filtered10_functional_stats.py", "--emapper-dir", f"{S}/func/em", "--ips-dir", f"{S}/func/ips",
         "--filtered", f"{S}/func/filtered", "--table", f"{S}/func/01_per_strain_functional_summary.tsv", "--orig", f"{S}/func/orig", "--out", ff],
         "SYNTHETIC: 5 proteins, p5 removed; IPS p1,p2(,p5); eggNOG p2,p3(,p5)", 0,
         "n=4 ips=2 egg=2 union=3 both=1 ips_only=1 egg_only=1 GO union=2 (InterPro 1, PANTHER 1) union_pct=75.0; SYN_B unchanged", check=ffcheck)

    def setops(_):
        rows = list(csv.DictReader(open(f"{I}/02_per_strain_functional_summary_final.tsv"), delimiter="\t"))
        bad = [r["strain"] for r in rows if not (int(r["union_hit"]) == int(r["both_tools"]) + int(r["ips_only"]) + int(r["egg_only"])
               and int(r["ips_hit"]) == int(r["both_tools"]) + int(r["ips_only"]) and int(r["egg_hit"]) == int(r["both_tools"]) + int(r["egg_only"]))]
        tot = sum(int(r["n_proteins"]) for r in rows)
        nf = sum(r["filtered"] == "yes" for r in rows)
        return (not bad and tot == 2618164 and len(rows) == 202 and nf == 10, f"rows={len(rows)} total={tot} filtered={nf} set-identity violations={len(bad)}")
    ok, det = setops(None)
    rec("A02", "real_small", "consistency of inputs/02_per_strain_functional_summary_final.tsv (union = both + ips_only + egg_only, ...)",
        "inputs/02_per_strain_functional_summary_final.tsv", "202 rows, 2,618,164 proteins, 10 filtered, 0 violations", "PASS" if ok else "FAIL", det, 0)

    # ---------------- 07 subphylum medians of union coverage -------------------
    import statistics as _st
    from decimal import Decimal as _D, ROUND_HALF_UP as _HU
    def _d1(x): return str(_D(repr(x)).quantize(_D("0.1"), rounding=_HU))
    rows_ = list(csv.DictReader(open(f"{I}/per_strain_release.tsv"), delimiter="\t"))
    med_ = {}
    for sp in ("Muc", "Mor", "Glo"):
        v = [100.0 * int(r["any"]) / int(r["n_release"]) for r in rows_ if r["subphylum"] == sp]
        med_[sp] = (len(v), _st.median(v))
    shown = {sp: _d1(m) for sp, (n, m) in med_.items()}
    ok = shown == {"Muc": "91.2", "Mor": "82.4", "Glo": "72.3"} and [med_[k][0] for k in ("Muc", "Mor", "Glo")] == [123, 57, 22]
    rec("A03", "real_small", "statistics.median of per-genome union coverage (100*any/n_release) by subphylum, one decimal (round half up)",
        "inputs/per_strain_release.tsv (counts)", "n 123/57/22; Muc 91.2, Mor 82.4, Glo 72.3 (manuscript display values)", "PASS" if ok else "FAIL",
        "; ".join(f"{sp} n={n} median={m:.6f} -> {shown[sp]}" for sp, (n, m) in med_.items()), 0)
    # SYNTHETIC: recheck_func_numbers.py writes unrounded medians; even n -> mean of the two middle values
    rs = os.path.join(W, "recheck_syn"); ips_d, em_d, pre_d, flt_d = (os.path.join(rs, x) for x in ("ips", "em", "pre", "flt"))
    for d_ in (ips_d, em_d, pre_d, flt_d): os.makedirs(d_)
    emcols = ["#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name",
              "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]
    spec = {"Glomus_SYN1": 700, "Glomus_SYN2": 722, "Glomus_SYN3": 723, "Glomus_SYN4": 750, "Mortierella_SYN1": 800, "Absidia_SYN1": 900}
    for lab, nh in spec.items():
        ids = [f"p{k:04d}" for k in range(1, 1001)]
        extra = [f"q{k:02d}" for k in range(1, 11)] if lab == "Absidia_SYN1" else []
        open(f"{pre_d}/{lab}_longest_isoform.faa", "w").write("".join(f">{x}\nMAAAA\n" for x in ids + extra))
        if extra: open(f"{flt_d}/{lab}_longest_isoform.filtered.faa", "w").write("".join(f">{x}\nMAAAA\n" for x in ids))
        lines = [f"{x}\tmd5\t100\t{'AntiFam' if (lab == 'Absidia_SYN1' and x in ('p0001', 'q01')) else 'Pfam'}\tPF1\td\t1\t50\t1e-5\tT\t01-01-2026\t-\t-\t-\t-\n"
                 for x in ids[:nh] + extra[:1]]
        open(f"{ips_d}/{lab}.interproscan.tsv", "w").write("".join(lines))
        open(f"{em_d}/{lab}.emapper.annotations", "w").write("## SYNTHETIC\n" + "\t".join(emcols) + "\n")
    def medcheck(_):
        t = {r["subphylum"]: r for r in csv.DictReader(open(f"{rs}/out/subphylum_medians.tsv"), delimiter="\t")}
        g = t["Glo"]
        ok = g["n_genomes"] == "4" and g["median_pct_any_unrounded"] == "72.250000" and g["median_pct_any_1dp"] == "72.3" and t["Mor"]["median_pct_any_1dp"] == "80.0"
        return ok, f"Glo n={g['n_genomes']} unrounded={g['median_pct_any_unrounded']} 1dp={g['median_pct_any_1dp']}; Mor {t['Mor']['median_pct_any_1dp']}; Muc {t['Muc']['median_pct_any_1dp']}"
    case("A04", "synthetic", [PY, "07_functional_annotation/recheck_func_numbers.py", "--ips-dir", ips_d, "--em-dir", em_d, "--prefilter-faa-dir", pre_d,
         "--filtered-dir", flt_d, "--out", f"{rs}/out"], "SYNTHETIC 6 genomes; Glo coverage 70.0/72.2/72.3/75.0", 0,
         "subphylum_medians.tsv: Glo n=4, unrounded 72.250000, one decimal 72.3 (half up); Mor 80.0", check=medcheck)

    # ---------------- 10 make_release / pack (synthetic release) ---------------
    rr = os.path.join(W, "rep_root"); fr = os.path.join(W, "func_release")
    for sub, src in (("01_202_longest_isoform_gtf", "gtf"), ("02_202_longest_isoform_proteins", "faa"), ("06_202_longest_isoform_gff3", "gff3")):
        os.makedirs(f"{rr}/{sub}")
        for f in glob.glob(f"{S}/release_check/ok/{src}/*"):
            shutil.copy(f, f"{rr}/{sub}/")
    os.makedirs(f"{fr}/SYN_A")
    open(f"{fr}/SYN_A/SYN_A_interproscan.tsv", "w").write("g1.t1\tSYNTHETIC\n")
    open(f"{fr}/SYN_A/SYN_A_emapper.annotations.tsv", "w").write("#SYNTHETIC\ng1.t1\n")
    rel = os.path.join(W, "release")
    case("P01", "synthetic", ["bash", "10_release_packaging/make_release_202.sh"], "SYNTHETIC 1-strain REP_ROOT/FUNC_RELEASE", 0,
         "release folder with 5 files, 0 ID mismatches, MD5SUMS_files.txt",
         env={"REP_ROOT": rr, "FUNC_RELEASE": fr, "OUT": rel},
         check=lambda o: (len(os.listdir(f"{rel}/SYN_A")) == 5 and os.path.exists(f"{rel}/MD5SUMS_files.txt")
                          and open(f"{rel}/00_release_check.tsv").read().splitlines()[1].split("\t")[4:] == ["0", "0", "0", "none"], "5 files; check row 0/0/0/none"))
    case("P02", "synthetic", ["bash", "10_release_packaging/make_release_202.sh"], "same, OUT exists", 1, "refuses",
         env={"REP_ROOT": rr, "FUNC_RELEASE": fr, "OUT": rel})
    case("P03", "synthetic", ["bash", "10_release_packaging/make_release_202.sh"], "SYNTHETIC 1 strain with STRICT=1", 3, "expectation 202/2,618,164 not met -> exit 3",
         env={"REP_ROOT": rr, "FUNC_RELEASE": fr, "OUT": rel + "_strict", "STRICT": "1"})
    smap = os.path.join(W, "smap.tsv"); open(smap, "w").write("strain\tsubphylum\nSYN_A\tMuc\n")
    case("P04", "synthetic", ["bash", "10_release_packaging/pack_release_202.sh"], "SYNTHETIC release (P01) + map SYN_A=Muc", 0,
         "Mucoromycotina/SYN_A.tar.gz with 5 members, gzip OK, MD5SUMS.txt and MD5SUMS_files.txt both kept",
         env={"RELEASE": rel, "SUBPHYLUM_MAP": smap},
         check=lambda o: (lambda t: (len([m for m in t.getmembers() if m.isfile()]) == 5 and os.path.exists(f"{rel}/MD5SUMS.txt") and os.path.exists(f"{rel}/MD5SUMS_files.txt"),
                                     ",".join(sorted(os.path.basename(m.name) for m in t.getmembers() if m.isfile()))))(tarfile.open(f"{rel}/Mucoromycotina/SYN_A.tar.gz")))
    case("P05", "synthetic", ["bash", "10_release_packaging/pack_release_202.sh"], "same, archives exist", 1, "refuses without OVERWRITE=1",
         env={"RELEASE": rel, "SUBPHYLUM_MAP": smap})
    case("P06", "synthetic", ["bash", "10_release_packaging/pack_release_202.sh"], "SYNTHETIC, STRICT=1 OVERWRITE=1", 3, "1 archive != 202 = 123/57/22 -> exit 3",
         env={"RELEASE": rel, "SUBPHYLUM_MAP": smap, "STRICT": "1", "OVERWRITE": "1"})
    case("P07", "synthetic", ["bash", "10_release_packaging/pack_release_202.sh"], "SYNTHETIC, map path does not exist", 2, "stops before packing",
         env={"RELEASE": rel, "SUBPHYLUM_MAP": os.path.join(W, "nope.tsv")})

    # ---------------- 08 published trees -------------------------------------
    if a.trees and have_dp:
        to = os.path.join(W, "trees_out")
        def tcheck(o):
            want = ["P_occ80_g90 vs G_occ80_g90: shared=192/199 RF=14 nRF=0.0352", "P_occ90_g90 vs G_occ90_g90: shared=194/199 RF=10 nRF=0.0251",
                    "P_occ80_g90 vs P_occ90_g90: shared=192/199 RF=14 nRF=0.0352", "G_occ80_g90 vs G_occ90_g90: shared=195/199 RF=8 nRF=0.0201"]
            got = [l.split(" strong")[0] for l in o.splitlines() if " vs " in l]
            return (got == want and os.path.exists(f"{to}/pairwise_summary.tsv"), "; ".join(got))
        case("T01", "real_small", [PY, "08_phylogenomics/compare_four_trees.py", "--tree-dir", a.trees, "--outdir", to],
             "four published treefiles (207 tips each)", 0, "202 ingroup tips after pruning 5 outgroups; P80vsG80 shared 192/199 nRF 0.0352 (values recorded from this run, not imposed)", check=tcheck)
        td = os.path.join(W, "trees_dup"); copy_plain(a.trees, td)
        t = open(f"{td}/G_occ90_g90.treefile").read()
        open(f"{td}/G_occ90_g90.treefile", "w").write(t.replace("Neurospora_crassa", "SYN_REMOVED_OUTGROUP", 1))
        case("T02", "synthetic", [PY, "08_phylogenomics/compare_four_trees.py", "--tree-dir", td, "--outdir", os.path.join(W, "trees_dup_out")],
             "copy of real trees, one outgroup label SYNTHETICALLY renamed", "nonzero", "stops (outgroup missing / tip sets differ)")
        tm = os.path.join(W, "trees_missing"); copy_plain(a.trees, tm); os.remove(f"{tm}/P_occ90_g90.treefile")
        case("T03", "synthetic", [PY, "08_phylogenomics/compare_four_trees.py", "--tree-dir", tm, "--outdir", os.path.join(W, "trees_missing_out")],
             "copy of real trees with one treefile removed", "nonzero", "stops (missing treefile)")
    else:
        why = "DendroPy not installed" if not have_dp else "no --trees given"
        for tid, cat, what in (("T01", "real_small", "compare_four_trees.py on the four treefiles"),
                               ("T02", "synthetic", "compare_four_trees.py with one outgroup label renamed (expect stop)"),
                               ("T03", "synthetic", "compare_four_trees.py with one treefile missing (expect stop)")):
            rec(tid, cat, what, "-", "see 08_phylogenomics/README.md", "NOT_RUN", why, "-")


    # ---------------- v2: R syntax, plotting scripts, naming/doc consistency --
    FIGD = "09_tables_and_figures/figures"
    rscript = a.rscript or shutil.which("Rscript")
    def py_has(exe, mods):
        if not exe: return False
        return subprocess.run([exe, "-c", "import " + ",".join(mods)], capture_output=True).returncode == 0
    plot_py = a.plot_python or next((x for x in (PY, shutil.which("python3")) if py_has(x, ["pandas", "matplotlib"])), None)
    rall = sorted(os.path.relpath(p, REPO) for p in glob.glob(os.path.join(REPO, "**", "*.[rR]"), recursive=True))
    rfiles, rhist = pub(rall), hist(rall)
    def r_bad(files):
        return [f for f in files if run([rscript, "-e", f"invisible(parse(file={f!r}, encoding='UTF-8'))"])[0]]
    if rscript:
        bad = r_bad(rfiles)
        rec("S04", "syntax", "Rscript -e parse(file=...) on every public (non-historical) *.R / *.r", f"{len(rfiles)} files",
            "all parse", "PASS" if not bad else "FAIL", f"{len(rfiles) - len(bad)}/{len(rfiles)} parse; failed: {bad}", len(bad))
    else:
        rec("S04", "syntax", "Rscript -e parse(file=...) on every public *.R / *.r", f"{len(rfiles)} files", "all parse", "NOT_RUN", "Rscript not found (use --rscript)", "-")
    if a.include_historical:
        if rscript:
            bad = r_bad(rhist)
            rec("S04h", "historical_audit", "Rscript -e parse(file=...) on historical/*.R / *.r", f"{len(rhist)} files",
                "all parse", "PASS" if not bad else "FAIL", f"{len(rhist) - len(bad)}/{len(rhist)} parse; failed: {bad}", len(bad))
        else:
            rec("S04h", "historical_audit", "Rscript parse on historical R files", f"{len(rhist)} files", "all parse", "NOT_RUN", "Rscript not found (use --rscript)", "-")
    plots_r = sorted(os.path.relpath(p, REPO) for p in glob.glob(os.path.join(REPO, FIGD, "*.R")))
    for i, f in enumerate(plots_r + [f"{FIGD}/plot_omark_contamination_filtering.py"], 1):
        tid = f"H{18 + i:02d}"
        exe = rscript if f.endswith(".R") else plot_py
        if not exe:
            rec(tid, "help", f"{'Rscript' if f.endswith('.R') else 'python3'} {f} --help", "-", "exit 0", "NOT_RUN",
                "Rscript not found" if f.endswith(".R") else "no python with pandas+matplotlib (use --plot-python)", "-"); continue
        case(tid, "help", [exe, f, "--help"], "-", 0, "usage printed, no plotting",
             check=lambda o: (o.strip().startswith("Usage:"), o.strip().splitlines()[0] if o.strip() else ""))

    # X01 relative calls ($here/..., $SD/...) and `source` targets exist
    import re as _re
    def current(rel):
        return not any(s in "/" + rel for s in ("/historical/", "/_versions/", "/evidence/", "/pre_refactor/")) and not rel.startswith(("inputs/", "tests/"))
    shs = [f for f in release_files("sh") if current(f)]
    miss = []
    for f in shs:
        for ln in open(os.path.join(REPO, f), encoding="utf-8", errors="replace"):
            if ln.lstrip().startswith("#"): continue
            for m in _re.finditer(r'\$\{?(?:here|SD)\}?/([A-Za-z0-9_./-]+\.(?:sh|py))', ln):
                if not os.path.exists(os.path.join(REPO, os.path.dirname(f), m.group(1))):
                    miss.append(f"{f}: {m.group(1)}")
    rec("X01", "consistency", "every $here/... and $SD/... script reference in current *.sh exists", f"{len(shs)} scripts",
        "0 missing", "PASS" if not miss else "FAIL", f"missing: {miss}" if miss else "0 missing", len(miss))

    # X02 repository paths in config examples and READMEs exist
    docs = sorted(set(glob.glob(os.path.join(REPO, "config", "*"))) | set(glob.glob(os.path.join(REPO, "**", "README.md"), recursive=True))
                  | set(glob.glob(os.path.join(REPO, "docs", "*.md"))))
    docs = [d for d in docs if current(os.path.relpath(d, REPO)) or os.path.relpath(d, REPO).startswith("tests/")]
    miss = []
    for d in docs:
        rel = os.path.relpath(d, REPO); base = os.path.dirname(d); in_map = False
        for ln in open(d, encoding="utf-8", errors="replace"):
            if ln.startswith("## Script names in run logs"): in_map = True
            if in_map or _re.search("historical", ln, _re.I): continue
            for m in _re.finditer(r'(?<![\w$/~.-])((?:\d\d_[a-z_]+|config|tests|docs|inputs)(?:/[A-Za-z0-9_.-]+)+\.(?:py|sh|R|tsv|example|md))(?![\w.])', ln):
                p = m.group(1)
                if not (os.path.exists(os.path.join(REPO, p)) or os.path.exists(os.path.join(base, p))):
                    miss.append(f"{rel}: {p}")
    rec("X02", "consistency", "repository paths named in config/*, README.md files and docs/*.md exist (lines about historical runs excluded)",
        f"{len(docs)} files", "0 missing", "PASS" if not miss else "FAIL", f"missing: {miss[:20]}" if miss else "0 missing", len(miss))

    # X03 no stale (pre-rename) script names in current files
    rmap = os.path.join(REPO, "docs", "RENAME_MAP.tsv")
    if os.path.exists(rmap):
        rows = list(csv.DictReader(open(rmap, encoding="utf-8"), delimiter="\t"))
        olds = {os.path.basename(r["Old path"]) for r in rows if r["Action"] == "rename"
                and os.path.basename(r["Old path"]) != os.path.basename(r["New path"]) and "/historical/" not in r["New path"]
                and not r["Old path"].startswith("uploaded:")}
        olddirs = {"04_structural_annotation/braker2_batch_v3/", "04_structural_annotation/braker3_batch_v4/",
                   "06_quality_and_filtering/omark_filtered_v6/", "06_quality_and_filtering/12_faa_12_gff_technical_validation/",
                   "ips_fix_20261011/fix_20261011_part2/"}
        pat = _re.compile(r'(?<![\w.-])(' + "|".join(_re.escape(o) for o in sorted(olds, key=len, reverse=True)) + r')(?![\w-])')
        hits = []; nf = 0
        for dp, dn, fn in os.walk(REPO):
            for fname in fn:
                rel = os.path.relpath(os.path.join(dp, fname), REPO)
                if not current(rel) or rel.startswith("docs/") and rel.endswith(".tsv"): continue
                if not fname.endswith((".sh", ".py", ".R", ".r", ".md", ".example", ".tsv", ".txt")): continue
                nf += 1; in_map = False
                for k, ln in enumerate(open(os.path.join(dp, fname), encoding="utf-8", errors="replace"), 1):
                    if ln.startswith("## Script names in run logs"): in_map = True
                    if in_map or _re.search(r"historical|\[release rename\]", ln, _re.I): continue
                    if pat.search(ln) or any(o in ln for o in olddirs):
                        hits.append(f"{rel}:{k}")
        rec("X03", "consistency", "old script names / old code-folder paths absent from current files (historical lines and rename tables excluded)",
            f"{nf} files, {len(olds)} old names", "0 hits", "PASS" if not hits else "FAIL", f"hits: {hits[:20]}" if hits else "0 hits", len(hits))
        bad = [r["New path"] for r in rows if r["Action"] in ("rename", "keep", "hold") and not r["New path"].startswith("-")
               and not os.path.exists(os.path.join(REPO, r["New path"]))]
        gone = [r["Old path"] for r in rows if r["Action"] == "rename" and not r["Old path"].startswith("uploaded:")
                and os.path.exists(os.path.join(REPO, r["Old path"]))]
        rec("X04", "consistency", "docs/RENAME_MAP.tsv: every new path exists; renamed old paths no longer exist", f"{len(rows)} rows",
            "0 problems", "PASS" if not (bad or gone) else "FAIL", f"missing new: {bad[:10]}; old still present: {gone[:10]}", len(bad) + len(gone))
    else:
        rec("X03", "consistency", "stale-name check", "-", "-", "FAIL", "docs/RENAME_MAP.tsv missing", "-")

    # X05 every current top-level script is named in its stage README
    unl = []
    for st in sorted(glob.glob(os.path.join(REPO, "[01][0-9]_*"))):
        readme = "".join(open(r, encoding="utf-8").read() for r in glob.glob(os.path.join(st, "**", "README.md"), recursive=True))
        for p in glob.glob(os.path.join(st, "**", "*"), recursive=True):
            rel = os.path.relpath(p, REPO)
            if os.path.isfile(p) and current(rel) and p.endswith((".py", ".sh", ".R")) and os.path.basename(p) not in readme:
                unl.append(rel)
    rec("X05", "consistency", "every current script (stages 01-10) is named in its stage README", "-", "0 unlisted",
        "PASS" if not unl else "FAIL", f"unlisted: {unl}" if unl else "0 unlisted", len(unl))

    # X06 figure index covers every plotting script; X07 required fields appear in the script
    idx = list(csv.DictReader(open(os.path.join(REPO, FIGD, "figure_index.tsv"), encoding="utf-8"), delimiter="\t"))
    srcs = list(csv.DictReader(open(os.path.join(REPO, FIGD, "plot_sources.tsv"), encoding="utf-8"), delimiter="\t"))
    listed = {r["script_path"] for r in idx} | {r["original_copy"] for r in srcs}
    sig = _re.compile(r"ggplot\(|ggsave\(|matplotlib|grDevices::pdf\(|savefig\(|<svg|svg_parts|\.svg\"")
    found = []
    for dp, dn, fn in os.walk(REPO):
        for fname in fn:
            rel = os.path.relpath(os.path.join(dp, fname), REPO)
            if rel.startswith(("tests/", "docs/")) or not fname.endswith((".py", ".R", ".r")): continue
            if sig.search(open(os.path.join(dp, fname), encoding="utf-8", errors="replace").read()): found.append(rel)
    unidx = [f for f in found if f not in listed]
    rec("X06", "figures", "every script with plotting calls is listed in figure_index.tsv or plot_sources.tsv", f"{len(found)} plotting scripts",
        "0 unlisted", "PASS" if not unidx else "FAIL", f"unlisted: {unidx}" if unidx else f"{len(found)} listed", len(unidx))
    nofield = []
    for r in idx:
        if not r["script_path"].startswith(FIGD + "/"): continue
        code = open(os.path.join(REPO, r["script_path"]), encoding="utf-8").read()
        for part in _re.split(r"[;:]", r["required_fields"]):
            for tok in _re.findall(r"[A-Za-z_\u4e00-\u9fff][\w\u4e00-\u9fff]*(?:\{[^}]*\}[\w]*)*", part):
                if tok in ("optional", "Fig2_tip_data", "tsv", "tip_metadata", "conflict_splits", "conflict_regions", "script", "checks",
                           "strains", "by", "subphylum", "n_buscos", "and") and tok not in code: continue
                t = tok.split("{")[0]
                if "{" in tok or tok.startswith(("omark_",)):
                    if "omark_" not in code: nofield.append(f"{r['script_path']}: {tok}")
                elif t and len(t) > 2 and t not in code and not t.isdigit():
                    nofield.append(f"{r['script_path']}: {t}")
    rec("X07", "figures", "required_fields in figure_index.tsv occur in the corresponding script", f"{sum(r['script_path'].startswith(FIGD) for r in idx)} scripts",
        "0 missing", "PASS" if not nofield else "FAIL", f"missing: {nofield[:20]}" if nofield else "all fields found", len(nofield))
    import hashlib as _h
    badsha = []
    for r in srcs:
        if r["original_copy"].startswith(FIGD):
            if _h.sha256(open(os.path.join(REPO, r["original_copy"]), "rb").read()).hexdigest() != r["uploaded_sha256"]: badsha.append(r["original_copy"])
            if _h.sha256(open(os.path.join(REPO, r["release_path"]), "rb").read()).hexdigest() != r["release_sha256"]: badsha.append(r["release_path"])
    rec("X08", "figures", "plot_sources.tsv: SHA256 of uploaded originals and release copies match the files", f"{len(srcs)} rows",
        "all match", "PASS" if not badsha else "FAIL", f"mismatch: {badsha}" if badsha else "all match", len(badsha))
    stmt = ("These scripts generate the initial plots used for figure preparation.\nThe publication figures were subsequently refined manually in\n"
            "Adobe Illustrator 2026. The scripts reproduce the initial plotting\nstage rather than the exact final appearance of the publication figures.")
    fr = open(os.path.join(REPO, FIGD, "README.md"), encoding="utf-8").read()
    forbidden = [w for w in ("pixel", "only fonts", "only the font", "unchanged data elements", ".ai files are included") if w in fr.lower()]
    rec("X09", "figures", "figures/README.md contains the Illustrator statement verbatim and no reproduction claims", "-", "present; 0 forbidden phrases",
        "PASS" if stmt in fr and not forbidden else "FAIL", f"statement {'present' if stmt in fr else 'MISSING'}; forbidden: {forbidden}", "-")

    # X10/X11 drafts that need no input (real run, written only into OUTDIR)
    for tid, scr, outs in (("X10", "plot_database_overlap_annotation_availability.R",
                            ["Fig1a_FINAL_proportional_circles_single_bar.pdf", "Fig1a_FINAL_proportional_circles_single_bar.tiff",
                             "Fig1a_FINAL_proportional_circles_single_bar_preview.png", "Fig1a_FINAL_counts.tsv"]),
                           ("X11", "plot_workflow_diagram.R", ["Fig_workflow_v4.pdf", "Fig_workflow_v4.png"])):
        od = os.path.join(O, "plots", scr[:-2])
        if not rscript:
            rec(tid, "figures", f"Rscript {FIGD}/{scr} --outdir=OUT", "none", "files written", "NOT_RUN", "Rscript not found", "-"); continue
        case(tid, "figures", [rscript, os.path.join(REPO, FIGD, scr), f"--outdir={od}"], "none (no input data)", 0,
             "draft files written to the test folder only", cwd=W,
             check=lambda o, od=od, outs=outs: (all(os.path.getsize(os.path.join(od, x)) > 0 for x in outs if os.path.exists(os.path.join(od, x)))
                                                and all(os.path.exists(os.path.join(od, x)) for x in outs), f"written: {sorted(os.listdir(od)) if os.path.isdir(od) else []}"))
    # X12-X14 SYNTHETIC inputs lacking one required column -> script stops before plotting
    syn = os.path.join(W, "plot_synthetic"); os.makedirs(syn)
    open(f"{syn}/SYN_busco_by_source_missing_col.tsv", "w").write(
        "subphylum\tspecies_no\tspecies_name_current\tgenome_BUSCO_C\tpublic_protein_BUSCO_C\nMucoromycotina\t1\tSYN_A\t95\t90\n")
    open(f"{syn}/SYN_busco12_missing_col.tsv", "w").write(
        "sample\tsample_order\tsubphylum\tmode\tlineage\tdataset_label\tmetric\tpercent\tcomplete_percent\nSYN_A\t1\tMucoromycotina\tgenome\tx\tG\tS\t90\t95\n")
    open(f"{syn}/SYN_omark_missing_col.tsv", "w").write("sample\tspecies_label\tsubphylum\nSYN_A\tSYN A\tMucoromycotina\n")
    for tid, exe, scr, inp, col in (("X12", rscript, "plot_busco_completeness_by_source.R", "SYN_busco_by_source_missing_col.tsv", "braker_protein_BUSCO_C"),
                                    ("X13", rscript, "plot_busco_genome_reference_study_12.R", "SYN_busco12_missing_col.tsv", "n_buscos"),
                                    ("X14", plot_py, "plot_omark_contamination_filtering.py", "SYN_omark_missing_col.tsv", "omark_before_C_pct")):
        if not exe:
            rec(tid, "synthetic", f"{scr} --input=SYNTHETIC", "SYNTHETIC", "nonzero", "NOT_RUN", "interpreter/packages missing", "-"); continue
        od = os.path.join(syn, "out_" + tid)
        case(tid, "synthetic", [exe, os.path.join(REPO, FIGD, scr), f"--input={syn}/{inp}", f"--outdir={od}"], f"SYNTHETIC {inp} (column {col} absent)", "nonzero",
             f"stops naming {col}; no figure written", cwd=W,
             check=lambda o, col=col, od=od: (col in o and not any(x.endswith((".pdf", ".png")) for x in (os.listdir(od) if os.path.isdir(od) else [])),
                                              o.strip().splitlines()[-1] if o.strip() else ""))
    # X15 caption constant in plot_protein_tree_annotation.R vs the published full trees
    if a.trees and have_dp:
        import dendropy
        from dendropy.calculate import treecompare
        tns = dendropy.TaxonNamespace()
        tp = dendropy.Tree.get(path=os.path.join(a.trees, "P_occ80_g90.treefile"), schema="newick", taxon_namespace=tns, preserve_underscores=True)
        tg = dendropy.Tree.get(path=os.path.join(a.trees, "G_occ80_g90.treefile"), schema="newick", taxon_namespace=tns, preserve_underscores=True)
        for t in (tp, tg): t.is_rooted = False; t.encode_bipartitions()
        def nsplit(t): return {b.split_bitmask if not b.split_bitmask & 1 else ~b.split_bitmask & ((1 << len(tns)) - 1)
                               for b in (e.bipartition for e in t.postorder_edge_iter()) if not b.is_trivial()}
        sp, sg = nsplit(tp), nsplit(tg); rf = treecompare.symmetric_difference(tp, tg)
        got = f"{len(sp & sg)}/{len(sp)}; nRF = {rf / (len(sp) + len(sg)):.3f}"
        note = "197/204; nRF = 0.034" in open(os.path.join(REPO, FIGD, "plot_protein_tree_annotation.R"), encoding="utf-8").read()
        rec("X15", "real_small", "dendropy: shared splits P_occ80_g90 vs G_occ80_g90, all 207 tips", "two published treefiles",
            "caption constant 197/204; nRF = 0.034", "PASS" if note and got == "197/204; nRF = 0.034" else "FAIL", got, 0)
    else:
        rec("X15", "real_small", "dendropy: shared splits P_occ80_g90 vs G_occ80_g90", "-", "197/204; nRF = 0.034", "NOT_RUN",
            "DendroPy not installed" if not have_dp else "no --trees given", "-")
    # X16 BRAKER2 best_by_compleasm scripts: historical names kept; only the config line differs from the pre-refactor copy
    hold = ["06_bbc.sh", "06_bbc_all.sh", "06b_bbc_merge.sh", "06c_bbc_full.sh", "09_rebuild_manifest.sh"]
    probs = []
    for h in hold:
        cur = os.path.join(REPO, "04_structural_annotation/braker2", h)
        old = os.path.join(REPO, "04_structural_annotation/historical/pre_refactor/braker2_batch_v3", h)
        if not (os.path.exists(cur) and os.path.exists(old)): probs.append(f"{h}: missing"); continue
        a1, b1 = open(old, encoding="utf-8").read().splitlines(), open(cur, encoding="utf-8").read().splitlines()
        import difflib as _dl
        for tag, i1, i2, j1, j2 in _dl.SequenceMatcher(None, a1, b1, autojunk=False).get_opcodes():
            if tag == "equal": continue
            for y in b1[j1:j2]:
                ok = (y.strip() == 'source "$here/config.sh"' or "[release refactor]" in y
                      or _re.search(r"\$(ENV_ROOT|BUSCO_LINEAGE_DIR|\{BBC_CWD)", y))
                if not ok: probs.append(f"{h}: unexpected line {y.strip()[:60]}")
    rec("X16", "consistency", "BRAKER2 best_by_compleasm scripts keep their historical names; vs the original only the config line and marked [release refactor] path lines differ", f"{len(hold)} scripts",
        "0 problems", "PASS" if not probs else "FAIL", f"{probs}" if probs else "only allowed lines differ", len(probs))
    # X17 public Markdown: no dead relative links or missing referenced files
    # (wording reminders for the authors are in delivery_notes/author_tools/check_public_wording.py, not part of this suite)
    mds = sorted(os.path.relpath(p, REPO) for p in glob.glob(os.path.join(REPO, "**", "*.md"), recursive=True))
    probs = []
    for rel in mds:
        txt = open(os.path.join(REPO, rel), encoding="utf-8").read()
        for m in _re.finditer(r"\]\(([^)#\s]+)(?:#[^)]*)?\)", txt):
            tgt = m.group(1)
            if tgt.startswith(("http://", "https://", "mailto:")): continue
            if not os.path.exists(os.path.normpath(os.path.join(REPO, os.path.dirname(rel), tgt))): probs.append(f"{rel}: dead link {tgt}")
    rec("X17", "consistency", "public *.md: no dead relative links", f"{len(mds)} files", "0 problems",
        "PASS" if not probs else "FAIL", f"{probs[:20]}" if probs else "0 problems", len(probs))
    # X18 config/phylogenomics_logged_runs.tsv vs the four IQ-TREE logs and the occ80 run script
    cfg = [r for r in csv.DictReader((l for l in open(os.path.join(REPO, "config/phylogenomics_logged_runs.tsv")) if not l.startswith("#")), delimiter="\t")]
    probs = []
    for r in cfg:
        lg = os.path.join(REPO, "08_phylogenomics/evidence/iqtree_logs", r["tree"] + ".log")
        if not os.path.exists(lg): probs.append(f"{r['tree']}: log missing"); continue
        txt = open(lg, encoding="utf-8", errors="replace").read()
        cmdl = next((l[len("Command:"):].split() for l in txt.splitlines() if l.startswith("Command:")), [])
        opt = lambda k: cmdl[cmdl.index(k) + 1] if k in cmdl else None
        sup = " ".join(x for x in ("-B " + str(opt("-B")) if opt("-B") else "", "-alrt " + str(opt("-alrt")) if opt("-alrt") else "",
                                    "-bnni" if "-bnni" in cmdl else "") if x)
        got = {"input_matrix": opt("-s"), "seq_type": opt("-st"), "model": opt("-m"), "support": sup, "threads": opt("-T"),
               "memory": opt("-mem"), "seed": opt("-seed"),
               "iqtree_version_in_log": (_re.search(r"IQ-TREE version (\S+)", txt) or [None, None])[1]}
        probs += [f"{r['tree']}: {k} config={r[k]} log={v}" for k, v in got.items() if r.get(k) != v]
    sh = open(os.path.join(REPO, "08_phylogenomics/run_iqtree_occ80.sh"), encoding="utf-8").read()
    if not _re.search(r'-m "\$MODEL" \\\n\s*-B 1000 -alrt 1000 -bnni \\\n\s*-T "\$THREADS" -mem "\$MEM" -seed "\$SEED"', sh): probs.append("run_iqtree_occ80.sh: tree command differs")
    if not all(x in sh for x in ("THREADS=20", "MEM=80G", "SEED=12345")): probs.append("run_iqtree_occ80.sh: occ80 threads/memory/seed differ")
    rec("X18", "consistency", "config/phylogenomics_logged_runs.tsv vs Command/version lines of 08_phylogenomics/evidence/iqtree_logs/*.log; run_iqtree_occ80.sh tree command",
        f"{len(cfg)} trees", "0 differences", "PASS" if not probs and len(cfg) == 4 else "FAIL", f"{probs}" if probs else "4 trees agree with their logs", len(probs))
    # P01-P07 plotting drafts that need inputs not in the repository
    PN = [("P01", "plot_busco_completeness_by_source.R", "--input=<per-record BUSCO table> --outdir=OUT", "per-record BUSCO table (TSV/CSV/XLSX)"),
          ("P02", "plot_umbelopsis_ramanniana_busco.R", "--input=Umbelopsis_ramanniana_BUSCO_completeness.tsv --outdir=OUT", "Umbelopsis_ramanniana_BUSCO_completeness*.tsv"),
          ("P03", "plot_genome_species_tree.R", "--tree=G_occ80_g90.treefile --metadata=tip_metadata.tsv --outdir=OUT", "tip_metadata.tsv (+ R packages ggtree, treeio)"),
          ("P04", "plot_omark_contamination_filtering.py", "--input=contamination_filtering_before_after_10strains_v2.tsv --outdir=OUT", "contamination_filtering_before_after_10strains_v2.tsv"),
          ("P05", "plot_busco_genome_vs_protein_by_subphylum.R", "--input=longest_isoform_busco_genome_protein_mucoromycota_202_v3-2.xlsx --outdir=OUT", "longest_isoform_busco_genome_protein_mucoromycota_202_v3-2.xlsx"),
          ("P06", "plot_busco_genome_reference_study_12.R", "--input=busco_12strains_GTR_long_v1.tsv --outdir=OUT", "busco_12strains_GTR_long_v1.tsv"),
          ("P07", "plot_protein_tree_annotation.R", "--input-dir=DIR --outdir=OUT", "DIR with P_*.treefile and Fig2_tip_data.tsv (+ ggtree, treeio, ggnewscale)")]
    for tid, scr, args_, need in PN:
        exe = "python3" if scr.endswith(".py") else "Rscript"
        rec(tid.replace("P", "D"), "plot_draft", f"{exe} {FIGD}/{scr} {args_}", need, "draft files written", "NOT_RUN", "input not in repository", "-")

    # ---------------- real-data regressions: documented, never executed here ---
    NR = [
        ("N01", "05 representative isoforms (213 BRAKER GTF)", "python3 05_representative_isoforms/filter_longest_isoforms_from_gtf.py --gtf-dir $BRAKER_GTF_DIR --protein-dir $BRAKER_PROTEIN_DIR --out-root $TMP/N01 && python3 05_representative_isoforms/repair_3_longest_isoform_samples.py --gtf-dir $BRAKER_GTF_DIR --protein-dir $BRAKER_PROTEIN_DIR --out-root $TMP/N01", "FAA/GTF identical (sha256) to historical 02_longest_isoform_proteins / 01_longest_isoform_gtf"),
        ("N02", "06 OMArk contig audit + filter on the prefilter view", "python3 06_quality_and_filtering/omark_filtering/audit_omark_contigs.py --omark-dir $OMARK_DIR --gtf-dir $PREFILTER_VIEW_GTF --genome-dir $GENOME_DIR --out $TMP/N02_v5 ; python3 06_quality_and_filtering/omark_filtering/filter_omark_contigs.py --v5-dirs $TMP/N02_v5 $V5_ACAU_DIR --gtf-dir $PREFILTER_VIEW_GTF --faa-dirs $PREFILTER_VIEW_FAA --genome-dir $GENOME_DIR --out $TMP/N02_v6", "10 filtered strains, 20,927 proteins removed, filtered FAA/GTF identical to omark_filtered_v6"),
        ("N03", "07 pre-filter functional summary", "python3 07_functional_annotation/summarize_functional_annotation.py --emapper-dir $FREEZE/emapper_v1_X --ips-dir $FREEZE/interproscan_v1_X --faa $PREFILTER_VIEW_FAA --out $TMP/N03", "01_per_strain_functional_summary.tsv identical to historical"),
        ("N04", "07 post-filter functional table", "python3 07_functional_annotation/recompute_filtered10_functional_stats.py --emapper-dir $FREEZE/emapper_v1_X --ips-dir $FREEZE/interproscan_v1_X --filtered $OMARK_FILTERED_FAA --table $TMP/N03/01_per_strain_functional_summary.tsv --orig $PREFILTER_VIEW_FAA --out $TMP/N04 && diff $TMP/N04/02_per_strain_functional_summary_final.tsv inputs/02_per_strain_functional_summary_final.tsv", "no diff"),
        ("N05", "07 manuscript recheck", "python3 07_functional_annotation/recheck_func_numbers.py --ips-dir $FUNC_RELEASE_IPS --em-dir $FUNC_RELEASE_EM --prefilter-faa-dir $PREFILTER_VIEW_FAA --filtered-dir $OMARK_FILTERED_FAA --out $TMP/N05 --subphylum-map $TMP/00_subphylum_map.tsv --fail-on-diff", "per_strain_release.tsv identical; subphylum_medians.tsv Glo 72.252888 (72.3 at one decimal); report.txt: median_any_Glo PASS against EXPECT 72.3; all items PASS"),
        ("N06", "05 GTF -> GFF3 (reconstructed wrapper)", "G2G=$AUGUSTUS_SCRIPTS/gtf2gff.pl GTF_DIR=$REP_ROOT/01_202_longest_isoform_gtf GFF3_DIR=$TMP/N06 bash 05_representative_isoforms/gtf_to_gff3.sh && diff -r $TMP/N06 $REP_ROOT/06_202_longest_isoform_gff3", "no diff apart from the two 00_* log files (historical gtf2gff.pl revision not recorded)"),
        ("N07", "10 release triplet check", "python3 10_release_packaging/check_release_inputs.py --labels inputs/release_labels_202.txt --faa-dir $REP_ROOT/02_202_longest_isoform_proteins --gtf-dir $REP_ROOT/01_202_longest_isoform_gtf --gff3-dir $REP_ROOT/06_202_longest_isoform_gff3 --expect-total 2618164 --out $TMP/N07.tsv", "202 OK, total 2,618,164"),
        ("N08", "10 release assembly + packaging", "REP_ROOT=$REP_ROOT FUNC_RELEASE=$FUNC_RELEASE OUT=$TMP/N08 STRICT=1 bash 10_release_packaging/make_release_202.sh && RELEASE=$TMP/N08 SUBPHYLUM_MAP=$TMP/00_subphylum_map.tsv STRICT=1 bash 10_release_packaging/pack_release_202.sh", "202 = 123/57/22 OK; per-file checksums equal to the Figshare checksum file (archive checksums may differ because of tar/gzip timestamps)"),
        ("N09", "08 matrix + IQ-TREE reruns", "see 08_phylogenomics/README.md (diagnose_and_select_markers.py -> align_trim_concat.py -> matrix_qc_and_column_filter.py -> run_iqtree_occ80.sh / compare_substitution_models.sh); needs the four supermatrices and locus lists (not bundled with this repository)", "four treefiles reproduce the published topologies"),
        ("N10", "04/06/07 heavy tools", "BRAKER2/BRAKER3, BUSCO, compleasm, OMArk, InterProScan, eggNOG-mapper batch runs", "not run by design (large databases, multi-day compute)"),
        ("N11", "06 12-genome public reference comparison", "REF12_BASE=... REF12_BUSCO_PROJECT=... BUSCO_DOWNLOAD_DIR=... bash 06_quality_and_filtering/reference12_validation/run_busco_ref12.sh ; REF12_BASE=... bash .../run_protein_comparison.sh", "tables equal to the historical 12-genome comparison"),
    ]
    for tid, what, cmd, exp in NR:
        rec(tid, "real_regression", cmd, what, exp, "NOT_RUN", "requires full data or heavy compute (not executed by this suite)", "-")

    # ---------------- environment + output ------------------------------------
    def ver(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True).stdout.strip().splitlines()[0]
        except Exception:
            return "n/a"
    env = {"python": sys.version.split()[0], "python_exe": PY, "platform": platform.platform(),
           "bash": ver(["bash", "--version"]), "perl": ver(["perl", "-e", 'print "perl $^V\n"']),
           "dendropy": __import__("dendropy").__version__ if have_dp else "not installed",
           "gtf2gff.pl": a.gtf2gff or "not given"}
    cols = ["id", "category", "command", "input", "expected", "actual", "exit_code", "status"]
    with open(os.path.join(O, "test_results.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in R:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
    cnt = {s: sum(r["status"] == s for r in R) for s in ("PASS", "FAIL", "NOT_RUN")}
    with open(os.path.join(O, "test_results.md"), "w") as fh:
        fh.write("# Test results\n\n" + "\n".join(f"- {k}: {v}" for k, v in env.items()) + "\n\n")
        fh.write(f"PASS {cnt['PASS']} / FAIL {cnt['FAIL']} / NOT_RUN {cnt['NOT_RUN']}\n\n")
        fh.write("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n")
        for r in R:
            fh.write("| " + " | ".join(str(r[c]).replace("|", "/") for c in cols) + " |\n")
    print(f"PASS {cnt['PASS']}  FAIL {cnt['FAIL']}  NOT_RUN {cnt['NOT_RUN']}  -> {O}/test_results.tsv")
    for r in R:
        if r["status"] != "PASS":
            print(f"  {r['status']}: {r['id']} {r['actual']}")
    sys.exit(1 if cnt["FAIL"] else 0)


if __name__ == "__main__":
    main()
