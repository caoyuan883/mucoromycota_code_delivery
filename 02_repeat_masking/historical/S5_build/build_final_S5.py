#!/usr/bin/env python3
"""合成表5 总表: 主表 + 软屏蔽核验, 并做名称归一、剔除非最终集、修正耗时列。"""
import os, re, csv, argparse
from collections import OrderedDict

KEY = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())

DROP = """Absidia_padenii Actinomortierella_ambigua Actinomortierella_wolfii
Amylomyces_rouxii Calcarisporiella_thermophila Cetraspora_pellucida
Dentiscutata_heterogama Endogone_sp._FLAS-F59071 Entrophospora_candida
Gamsiella_multidivaricata Haplosporangium_bisporale Jimgerdemannia_flammicorona
Jimgerdemannia_lactiflua Linnemannia_camargensis Mortierella_elongata
Mortierella_gamsii Mortierella_globulifera Mortierella_humilis Mucor_sp._PG5414
Podila_humilis Podila_minutissima Racocetra_fulgida Racocetra_persica
Rhizophagus_cerebriforme Rhizophagus_diaphanus Rhizophagus_diaphanu
Rhizophagus_sp._LPA8_3_AMF Rhizophagus_sp._LPA8_4_AMF
Rhizophagus_sp._MUCL_43196 Rhizopus_sp._PT2906 Scutellospora_calospora""".split()
DROPK = {KEY(x) for x in DROP}
RENAME = {KEY("Benjaminiella_poitrasii_JGI"): "Benjaminiella_poitrasii"}

DROP_COLS = ["satellite_n", "satellite_bp", "satellite_pct", "repeatmodeler_version",
             "rmod_search_engine", "rmod_dependencies", "repeatmasker_version",
             "softmask_xsmall", "pipeline_wallclock", "pipeline_wallclock_h", "lib_arg"]

def hms2h(x):
    try:
        h, m, s = map(int, str(x).split(":")); return h + m / 60 + s / 3600
    except Exception:
        return None

a = argparse.ArgumentParser()
a.add_argument("--master", default="S5_output/S5_final_table_v2.tsv")
a.add_argument("--softmask", default="softmask_check2.tsv")
a.add_argument("--out", default="S5_output/TableS5_final.tsv")
A = a.parse_args()

rows = list(csv.DictReader(open(A.master), delimiter="\t"))
print("主表行数:", len(rows))

sm_name, sm_len = {}, {}
for ln in open(A.softmask):
    p = ln.rstrip("\n").split("\t")
    if len(p) < 5 or not p[1].isdigit():
        continue
    rec = dict(lowercase_bp=int(p[1]), uppercase_bp=int(p[2]),
               N_total=int(p[3]), seq_bp=int(p[4]), file_name=p[0])
    sm_name[KEY(p[0])] = rec
    sm_len.setdefault(int(p[4]), rec)
print("软屏蔽记录:", len(sm_name))

out, stat = [], OrderedDict()
def bump(k): stat[k] = stat.get(k, 0) + 1

for r in rows:
    kk = KEY(r["sample"])
    if kk in DROPK:
        bump("剔除_非最终集"); continue
    if kk in RENAME:
        r["sample"] = RENAME[kk]; kk = KEY(r["sample"])

    tl = int(float(r.get("total_length_bp") or 0))
    tm = int(float(r.get("total_masked_bp") or 0))

    rec = sm_name.get(kk)
    how = "name"
    if rec is None:
        rec = next((v for k2, v in sm_name.items()
                    if k2.startswith(kk) or kk.startswith(k2)), None)
        how = "name_prefix" if rec else how
    if rec is None and tl:
        rec = sm_len.get(tl); how = "length" if rec else how

    if rec:
        lc, sb = rec["lowercase_bp"], rec["seq_bp"]
        r["softmask_file"] = rec["file_name"]
        r["softmask_match_by"] = how
        r["softmask_seq_bp"] = sb
        r["softmask_lowercase_bp"] = lc
        r["softmask_lowercase_pct"] = round(lc / sb * 100, 3) if sb else ""
        r["softmask_N_total"] = rec["N_total"]
        r["length_check"] = "OK" if sb == tl else "MISMATCH"
        r["softmask_vs_tbl_diff_bp"] = lc - tm
        r["softmask_vs_tbl_diff_pct"] = round((lc - tm) / tm * 100, 3) if tm else ""
        r["softmask_status"] = ("OK" if tm and abs(lc - tm) / tm < 0.01
                                else ("HIGHER" if lc > tm else "LOWER"))
        bump("软屏蔽_" + r["softmask_status"]); bump("长度_" + r["length_check"])
    else:
        for c in ("softmask_file", "softmask_match_by", "softmask_seq_bp",
                  "softmask_lowercase_bp", "softmask_lowercase_pct", "softmask_N_total",
                  "length_check", "softmask_vs_tbl_diff_bp",
                  "softmask_vs_tbl_diff_pct", "softmask_status"):
            r[c] = ""
        r["softmask_status"] = "FILE_MISSING"; bump("软屏蔽_FILE_MISSING")

    h = hms2h(r.get("rmod_time"))
    th = r.get("rmod_threads")
    r["repeatmodeler_time_h"] = round(h, 3) if h else ""
    try:
        r["cpu_core_hours"] = round(h * float(th), 1) if h and th else ""
    except Exception:
        r["cpu_core_hours"] = ""
    for c in DROP_COLS:
        r.pop(c, None)
    out.append(r); bump("保留")

cols = []
for r in out:
    for c in r:
        if c not in cols: cols.append(c)
with open(A.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in out: w.writerow(r)

print("\n===== 汇总 =====")
for k2, v in stat.items(): print(f"  {k2:24s} {v}")
tot = sum(float(r["cpu_core_hours"]) for r in out if r["cpu_core_hours"] != "")
print(f"  总算力(核·时)            {tot:.1f}")
print(f"[写出] {A.out}  ({len(out)} 行 x {len(cols)} 列)")
