#!/usr/bin/env python3
"""
从原始结果重新计算正文"功能注释的覆盖度"一节中的全部数字，并与正文逐项对照。
只读，不修改任何输入文件。

用法（release）：python3 recheck_func_numbers.py --help
（历史用法：cd ~/Desktop/interproscan_eggnog && python3 recheck_func_numbers.py）

输入（默认路径，与 recompute_filtered10_functional_stats.py 相同）：
  InterProScan : ~/Desktop/interproscan_eggnog/func_freeze_2/interproscan_v1_20261011/<株名>.interproscan.tsv.gz
  eggNOG-mapper: ~/Desktop/interproscan_eggnog/func_freeze_2/emapper_v1_20261011/<株名>.emapper.annotations(.gz)
  原始蛋白     : ~/Desktop/Mucoromycota_BRAKER_longest_isoform/02_202_longest_isoform_proteins/<株名>_longest_isoform.faa
  过滤版蛋白   : ~/Desktop/10_omark_filtered_faa/<株名>_longest_isoform.filtered.faa
输出：
  ~/Desktop/interproscan_eggnog/recheck_func_numbers/per_strain_release.tsv
  ~/Desktop/interproscan_eggnog/recheck_func_numbers/filtered10_before_after.tsv
  ~/Desktop/interproscan_eggnog/recheck_func_numbers/report.txt   （与正文逐项对照，PASS / DIFF）

口径：
  InterProScan 注释 = 该蛋白在 TSV 中至少有一行（任一成员库）
  eggNOG-mapper 注释 = 该蛋白出现在 .emapper.annotations 中
  GO：InterProScan 第 14 列非空且不为 "-"；eggNOG-mapper 的 GOs 列不为 "-"
  KEGG KO：eggNOG-mapper 的 KEGG_ko 列不为 "-"
  AntiFam：InterProScan 第 4 列为 AntiFam 的蛋白
  发布版：10 株用过滤版 FAA 中的蛋白 ID，其余 192 株用原始 FAA
"""
import os, sys, gzip, glob, statistics, argparse, csv

# [release refactor] The historical defaults above (~/Desktop/...) and the IPS_DIR/EM_DIR/
# FAA_DIR/FILT_DIR/OUT_DIR environment overrides are replaced by explicit arguments.
# IMPORTANT: --prefilter-faa-dir must contain the PRE-filter proteins of all 202 genomes
# (192 unchanged + 10 pre-filter originals). The directory 02_202_longest_isoform_proteins
# named in the historical docstring now holds the post-filter files for the 10 genomes;
# use 05_representative_isoforms/build_prefilter_view.py to assemble a pre-filter view.
# All calculations are unchanged (median = statistics.median). One EXPECT value changed:
# median_any_Glo 72.2 -> 72.3 (see the comment in EXPECT; historical original in historical/pre_refactor/).
_ap = argparse.ArgumentParser(description="Recompute functional-coverage numbers and compare with manuscript values (read-only).")
_ap.add_argument("--ips-dir", required=True, help="frozen InterProScan TSV directory (historical: func_freeze_2/interproscan_v1_20261011)")
_ap.add_argument("--em-dir", required=True, help="frozen eggNOG directory (historical: func_freeze_2/emapper_v1_20261011)")
_ap.add_argument("--prefilter-faa-dir", required=True, help="PRE-filter <strain>_longest_isoform.faa for all 202 genomes")
_ap.add_argument("--filtered-dir", required=True, help="<strain>_longest_isoform.filtered.faa for the post-filtered genomes")
_ap.add_argument("--out", required=True, help="output directory (must not exist or be empty)")
_ap.add_argument("--subphylum-map", default=None,
                 help="optional strain<TAB>subphylum(Muc/Mor/Glo) table; reports disagreements with the built-in genus lists (exit 4)")
_ap.add_argument("--fail-on-diff", action="store_true", help="exit 5 if any item is DIFF (default: report only, as historical)")
_args = _ap.parse_args()
IPS_DIR, EM_DIR, FAA_DIR, FILT_DIR, OUT = _args.ips_dir, _args.em_dir, _args.prefilter_faa_dir, _args.filtered_dir, _args.out
for _d in (IPS_DIR, EM_DIR, FAA_DIR, FILT_DIR):
    if not os.path.isdir(_d): sys.exit(f"ERROR: directory not found: {_d}")
if os.path.isdir(OUT) and os.listdir(OUT): sys.exit(f"ERROR: output directory not empty: {OUT}")

GLO = set("Acaulospora Ambispora Appendiculispora Dentiscutata Diversispora Entrophospora "
          "Funneliformis Geosiphon Gigaspora Glomus Oehlia Paraglomus Rhizophagus".split())
MOR = set("Actinomortierella Cryptocalciella Dissophora Entomortierella Gamsiella Gryganskiella "
          "Haplosporangium Linnemannia Lobosporangium Lunasporangiospora Modicella Mortierella Podila".split())

# 正文中的数字（百分数按一位小数比较）
EXPECT = {
    "release_total": 2618164, "prefilter_total": 2639091, "removed_total": 20927,
    "n_filtered_strains": 10,
    "pct_ips": 76.4, "pct_em": 81.9, "pct_any": 85.0, "pct_go_any": 65.1, "pct_ko": 38.2,
    "pct_go_em": 36.7, "pct_go_ips": 61.9,
    # [release refactor] median_any_Glo was 72.2; set to the manuscript display value 72.3.
    # The recomputed value (statistics.median = 72.252888) is still compared at one decimal;
    # evidence/recheck_report.txt was written with the earlier value 72.2 and is kept unchanged.
    "median_any_Muc": 91.2, "median_any_Mor": 82.4, "median_any_Glo": 72.3,
    "n_Muc": 123, "n_Mor": 57, "n_Glo": 22,
    "lowest6_all_Glo": True, "lowest6_min": 60.9, "lowest6_max": 67.4,
    "antifam_kept": 66, "antifam_kept_strains": 38, "antifam_removed": 10, "antifam_ratio": 19,
    "Rhizophagus_prolifer_before": 76.7, "Rhizophagus_prolifer_after": 73.7,
    "Acaulospora_colombiana_removed_em": 54.9,
    "Acaulospora_colombiana_before": 67.1, "Acaulospora_colombiana_after": 74.3,
}


def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def faa_ids(p):
    ids = []
    with opn(p) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def first(patterns):
    for pat in patterns:
        hit = glob.glob(pat)
        if hit:
            return hit[0]
    return None


def read_ips(p):
    hit, go, antifam = set(), set(), set()
    with opn(p) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            pid = c[0]
            hit.add(pid)
            if c[3] == "AntiFam":
                antifam.add(pid)
            if len(c) > 13 and c[13] not in ("", "-"):
                go.add(pid)
    return hit, go, antifam


def read_em(p):
    hit, go, ko = set(), set(), set()
    cols = None
    with opn(p) as f:
        for line in f:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if line.startswith("#query"):
                cols = {k: i for i, k in enumerate(c)}
                continue
            if cols is None or not c[0]:
                continue
            pid = c[0]
            hit.add(pid)
            if c[cols["GOs"]] not in ("", "-"):
                go.add(pid)
            if c[cols["KEGG_ko"]] not in ("", "-"):
                ko.add(pid)
    if cols is None:
        sys.exit(f"ERROR: eggNOG 表头 #query 未找到：{p}")
    return hit, go, ko


def subphylum(s):
    g = s.split("_")[0]
    return "Glo" if g in GLO else "Mor" if g in MOR else "Muc"


def pct(a, b):
    return 100.0 * a / b if b else float("nan")


os.makedirs(OUT, exist_ok=True)
strains = sorted(os.path.basename(p)[:-len("_longest_isoform.faa")]
                 for p in glob.glob(f"{FAA_DIR}/*_longest_isoform.faa"))
if len(strains) != 202:
    print(f"WARNING: 原始 FAA 数为 {len(strains)}，不是 202")

rows, filt_rows = [], []
T = dict(pre=0, rel=0, ips=0, em=0, any=0, go_any=0, go_ips=0, go_em=0, ko=0,
         af_kept=0, af_removed=0, removed=0)
af_strains = set()
for s in strains:
    orig = faa_ids(f"{FAA_DIR}/{s}_longest_isoform.faa")
    ips_p = first([f"{IPS_DIR}/{s}.interproscan.tsv.gz", f"{IPS_DIR}/{s}.interproscan.tsv"])
    em_p = first([f"{EM_DIR}/{s}.emapper.annotations", f"{EM_DIR}/{s}.emapper.annotations.gz"])
    if not ips_p or not em_p:
        sys.exit(f"ERROR: 缺少结果文件：{s}  ips={ips_p}  emapper={em_p}")
    ih, ig, ia = read_ips(ips_p)
    eh, eg, ek = read_em(em_p)
    fp = f"{FILT_DIR}/{s}_longest_isoform.filtered.faa"
    kept = faa_ids(fp) if os.path.exists(fp) else orig
    so, sk = set(orig), set(kept)
    if not sk <= so:
        sys.exit(f"ERROR: {s} 的过滤版含原始集中没有的 ID")
    if len(sk) != len(kept):
        sys.exit(f"ERROR: {s} 的过滤版有重复 ID")
    removed = so - sk
    n = len(sk)
    any_ = (ih | eh) & sk
    r = dict(strain=s, subphylum=subphylum(s), filtered=os.path.exists(fp), n_prefilter=len(so),
             n_release=n, n_removed=len(removed),
             ips=len(ih & sk), em=len(eh & sk), any=len(any_),
             go_any=len((ig | eg) & sk), go_ips=len(ig & sk), go_em=len(eg & sk), ko=len(ek & sk),
             antifam=len(ia & sk), antifam_removed=len(ia & removed))
    r["pct_any"] = pct(r["any"], n)
    r["pct_ips"] = pct(r["ips"], n)
    r["pct_em"] = pct(r["em"], n)
    rows.append(r)
    T["pre"] += len(so); T["rel"] += n; T["removed"] += len(removed)
    for k in ("ips", "em", "any", "go_any", "go_ips", "go_em", "ko"):
        T[k] += r[k]
    T["af_kept"] += r["antifam"]; T["af_removed"] += r["antifam_removed"]
    if r["antifam"]:
        af_strains.add(s)
    if r["filtered"]:
        before_any = len((ih | eh) & so)
        filt_rows.append(dict(strain=s, n_before=len(so), n_after=n, n_removed=len(removed),
                              pct_any_before=pct(before_any, len(so)), pct_any_after=r["pct_any"],
                              pct_em_removed=pct(len(eh & removed), len(removed)),
                              pct_em_kept=pct(len(eh & sk), n),
                              antifam_before=len(ia & so), antifam_after=r["antifam"]))

with open(f"{OUT}/per_strain_release.tsv", "w") as f:
    keys = list(rows[0].keys())
    f.write("\t".join(keys) + "\n")
    for r in rows:
        f.write("\t".join(f"{r[k]:.2f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
with open(f"{OUT}/filtered10_before_after.tsv", "w") as f:
    keys = list(filt_rows[0].keys()) if filt_rows else []
    f.write("\t".join(keys) + "\n")
    for r in sorted(filt_rows, key=lambda x: x["strain"]):
        f.write("\t".join(f"{r[k]:.2f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")

GOT = {
    "release_total": T["rel"], "prefilter_total": T["pre"], "removed_total": T["removed"],
    "n_filtered_strains": len(filt_rows),
    "pct_ips": pct(T["ips"], T["rel"]), "pct_em": pct(T["em"], T["rel"]),
    "pct_any": pct(T["any"], T["rel"]), "pct_go_any": pct(T["go_any"], T["rel"]),
    "pct_ko": pct(T["ko"], T["rel"]), "pct_go_em": pct(T["go_em"], T["rel"]),
    "pct_go_ips": pct(T["go_ips"], T["rel"]),
}
for sp in ("Muc", "Mor", "Glo"):
    v = [r["pct_any"] for r in rows if r["subphylum"] == sp]
    GOT[f"n_{sp}"] = len(v)
    GOT[f"median_any_{sp}"] = statistics.median(v) if v else float("nan")
low6 = sorted(rows, key=lambda r: r["pct_any"])[:6]
GOT["lowest6_all_Glo"] = all(r["subphylum"] == "Glo" for r in low6)
GOT["lowest6_min"] = low6[0]["pct_any"]; GOT["lowest6_max"] = low6[-1]["pct_any"]
GOT["antifam_kept"] = T["af_kept"]; GOT["antifam_kept_strains"] = len(af_strains)
GOT["antifam_removed"] = T["af_removed"]
GOT["antifam_ratio"] = (T["af_removed"] / T["removed"]) / (T["af_kept"] / T["rel"]) if T["af_kept"] and T["removed"] else float("nan")
fr = {r["strain"]: r for r in filt_rows}
for s in ("Rhizophagus_prolifer", "Acaulospora_colombiana"):
    if s in fr:
        GOT[f"{s}_before"] = fr[s]["pct_any_before"]; GOT[f"{s}_after"] = fr[s]["pct_any_after"]
GOT["Acaulospora_colombiana_removed_em"] = fr.get("Acaulospora_colombiana", {}).get("pct_em_removed", float("nan"))


def same(k, e, g):
    if isinstance(e, bool):
        return e == g
    if isinstance(e, int) and k != "antifam_ratio":
        return e == g
    if k == "antifam_ratio":
        return round(g) == e
    return round(g + 1e-9, 1) == e


lines = ["项目\t正文\t复算\t结果"]
bad = 0
for k, e in EXPECT.items():
    g = GOT.get(k, "NA")
    ok = g != "NA" and same(k, e, g)
    bad += (not ok)
    gs = f"{g:.2f}" if isinstance(g, float) else str(g)
    lines.append(f"{k}\t{e}\t{gs}\t{'PASS' if ok else 'DIFF'}")
lines.append("")
lines.append("最低 6 株：" + "; ".join(f"{r['strain']} {r['pct_any']:.2f}% ({r['subphylum']})" for r in low6))
lines.append("10 株被移除蛋白的 eggNOG 注释率（%）：" + "; ".join(
    f"{r['strain']} {r['pct_em_removed']:.1f}" for r in sorted(filt_rows, key=lambda x: x["pct_em_removed"])))
lines.append(f"\n共 {len(EXPECT)} 项，DIFF {bad} 项")
rep = "\n".join(lines)
open(f"{OUT}/report.txt", "w").write(rep + "\n")

# [release refactor] unrounded subphylum medians (statistics.median of per-genome pct_any;
# even n -> mean of the two middle values) and their one-decimal display value (round half up).
from decimal import Decimal, ROUND_HALF_UP
with open(f"{OUT}/subphylum_medians.tsv", "w") as _fh:
    _fh.write("subphylum\tn_genomes\tmedian_pct_any_unrounded\tmedian_pct_any_1dp\n")
    for _sp in ("Muc", "Mor", "Glo"):
        _m = GOT[f"median_any_{_sp}"]
        _d = Decimal(repr(_m)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP) if _m == _m else "NA"
        _fh.write(f"{_sp}\t{GOT['n_' + _sp]}\t{_m:.6f}\t{_d}\n")
print(rep)
print(f"\n输出目录：{OUT}")

# [release refactor] optional cross-check of the genus-list subphylum assignment
if _args.subphylum_map:
    with open(_args.subphylum_map, encoding="utf-8") as _fh:
        _m = {r["strain"]: r["subphylum"] for r in csv.DictReader(_fh, delimiter="\t")}
    _dis = [(r["strain"], r["subphylum"], _m.get(r["strain"], "MISSING")) for r in rows if _m.get(r["strain"]) != r["subphylum"]]
    with open(f"{OUT}/subphylum_crosscheck.tsv", "w") as _fh:
        _fh.write("strain\tgenus_list\tmap\n")
        for _r in _dis: _fh.write("\t".join(_r) + "\n")
    print(f"subphylum cross-check: {len(_dis)} disagreement(s)")
    if _dis: sys.exit(4)
if _args.fail_on_diff and bad: sys.exit(5)
