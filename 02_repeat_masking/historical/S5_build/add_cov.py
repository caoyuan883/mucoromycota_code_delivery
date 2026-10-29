#!/usr/bin/env python3
# 补回 RepeatModeler 采样覆盖度列
import os, re, csv, argparse

def key(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())

a = argparse.ArgumentParser()
a.add_argument("--table", default="S5_output/S5_final_table.tsv")
a.add_argument("--rmod-dir", required=True)
a.add_argument("--out", default="S5_output/S5_final_table_v2.tsv")
A = a.parse_args()

logs = {}
d = os.path.expanduser(A.rmod_dir)
for fn in os.listdir(d):
    if fn.endswith(".log"):
        logs[key(re.sub(r"[-_]?rmod\.log$", "", fn))] = os.path.join(d, fn)

rows = list(csv.DictReader(open(A.table), delimiter="\t"))
hit = 0
for r in rows:
    k = key(r["sample"])
    p = logs.get(k) or next((v for kk, v in logs.items()
                             if kk.startswith(k) or k.startswith(kk)), None)
    cov_bp = cov_pct = ""
    if p:
        t = open(p, "r", errors="replace").read()
        m = re.findall(r"Input Database Coverage:\s*(\d+)\s*bp out of\s*(\d+)\s*bp\s*\(\s*([\d.]+)\s*%", t)
        if m:
            cov_bp, _, cov_pct = m[-1]
            hit += 1
    r["sampled_bp"] = cov_bp
    r["sampling_coverage_pct"] = cov_pct

cols = list(rows[0].keys())
with open(A.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows: w.writerow(r)
print(f"抓到覆盖度的样本: {hit}/{len(rows)}")
print(f"[写出] {A.out} ({len(rows)} x {len(cols)})")
