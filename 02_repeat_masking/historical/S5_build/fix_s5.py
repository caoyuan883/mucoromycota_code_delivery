#!/usr/bin/env python3
# S5 后处理: 修正 wallclock、剔除非最终数据集的菌株、删除全空列
import os, re, csv, argparse
from datetime import datetime, timedelta

DROP = """Absidia_padenii Actinomortierella_ambigua Actinomortierella_wolfii
Amylomyces_rouxii Calcarisporiella_thermophila Cetraspora_pellucida
Dentiscutata_heterogama Endogone_sp._FLAS-F59071 Entrophospora_candida
Gamsiella_multidivaricata Haplosporangium_bisporale Jimgerdemannia_flammicorona
Jimgerdemannia_lactiflua Linnemannia_camargensis Mortierella_elongata
Mortierella_gamsii Mortierella_globulifera Mucor_sp._PG5414 Podila_humilis
Podila_minutissima Racocetra_fulgida Racocetra_persica Rhizophagus_cerebriforme
Rhizophagus_diaphanus Rhizophagus_sp._LPA8_3_AMF Rhizophagus_sp._LPA8_4_AMF
Rhizophagus_sp._MUCL_43196 Rhizopus_sp._PT2906 Scutellospora_calospora""".split()

def key(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())

DROPK = {key(x) for x in DROP}

def wall(path):
    ts = re.findall(r"\b(\d{1,2}:\d{2}:\d{2})\b", open(path, "r", errors="replace").read())
    if len(ts) < 2:
        return None
    f = datetime.strptime(ts[0], "%H:%M:%S")
    l = datetime.strptime(ts[-1], "%H:%M:%S")
    d = l - f
    if d.total_seconds() < 0:
        d += timedelta(days=1)
    h, r = divmod(int(d.total_seconds()), 3600)
    m, s = divmod(r, 60)
    return f"{h:02d}:{m:02d}:{s:02d}", round(d.total_seconds() / 3600, 3)

a = argparse.ArgumentParser()
a.add_argument("--master", default="S5_output/S5_master_table.tsv")
a.add_argument("--log-dir")
a.add_argument("--out", default="S5_output/S5_final_table.tsv")
A = a.parse_args()

rows = list(csv.DictReader(open(A.master), delimiter="\t"))
print("输入行数:", len(rows))

logs = {}
if A.log_dir:
    d = os.path.expanduser(A.log_dir)
    for fn in os.listdir(d):
        if fn.endswith((".log", ".txt")):
            logs[key(re.sub(r"\.(log|txt)$", "", fn))] = os.path.join(d, fn)

kept, dropped = [], []
for r in rows:
    k = key(r["sample"])
    if k in DROPK:
        dropped.append(r["sample"]); continue
    if logs:
        p = logs.get(k) or next((v for kk, v in logs.items()
                                 if kk.startswith(k) or k.startswith(kk)), None)
        if p:
            w = wall(p)
            if w:
                r["pipeline_wallclock"], r["pipeline_wallclock_h"] = w
    kept.append(r)

print("剔除:", len(dropped), "| 保留:", len(kept))
for s in dropped: print("   -", s)

cols = [c for c in kept[0].keys() if any(r.get(c) not in (None, "") for r in kept)]
empty = [c for c in kept[0].keys() if c not in cols]
if empty: print("删除全空列:", ", ".join(empty))

with open(A.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in kept: w.writerow(r)
print(f"[写出] {A.out}  ({len(kept)} x {len(cols)})")
