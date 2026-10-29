#!/usr/bin/env python3
"""校验软屏蔽产物: 统计各基因组 FASTA 中的小写碱基数, 与表5 的 total_masked_bp 比对。"""
import os, re, csv, argparse

KEY = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())
LOW = [b"a", b"c", b"g", b"t"]

def strip_name(fn):
    b = os.path.basename(fn)
    for suf in (".masked", ".gz"):
        if b.endswith(suf): b = b[: -len(suf)]
    b = re.sub(r"\.(fna|fa|fasta)$", "", b)
    b = re.sub(r"_genomic$", "", b)
    return re.sub(r"[._-]+$", "", b) 

def scan(path):
    seq = low = upN = lowN = 0    
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                continue
            line = line.strip()
            if not line:
                continue
            seq += len(line)
            for c in LOW:
                low += line.count(c)
            upN += line.count(b"N")
            lowN += line.count(b"n")
    return seq, low, upN, lowN

a = argparse.ArgumentParser()
a.add_argument("--genome-dir", required=True)
a.add_argument("--table", required=True)
a.add_argument("--out", default="softmask_validation.tsv")
A = a.parse_args()

gdir = os.path.expanduser(A.genome_dir)
files = {}
for fn in sorted(os.listdir(gdir)):
    p = os.path.join(gdir, fn)
    if os.path.isfile(p) and not fn.startswith("."):
        files.setdefault(KEY(strip_name(fn)), p)
print(f"发现基因组文件 {len(files)} 份")

rows = list(csv.DictReader(open(os.path.expanduser(A.table)), delimiter="\t"))
print(f"表5 记录 {len(rows)} 行")

out, tally = [], {}
for i, r in enumerate(rows, 1):
    s = r["sample"]; k = KEY(s)
    p = files.get(k) or next((v for kk, v in files.items()
                              if kk.startswith(k) or k.startswith(kk)), None)
    try:
        tm = int(float(r.get("total_masked_bp") or 0))
    except ValueError:
        tm = 0
    if not p:
        out.append(dict(sample=s, file="", seq_bp="", lowercase_bp="", upper_N="",
                        lower_n="", table_masked_bp=tm, diff_bp="", diff_pct="",
                        status="FILE_MISSING"))
        tally["FILE_MISSING"] = tally.get("FILE_MISSING", 0) + 1
        continue
    seq, low, upN, lowN = scan(p)
    diff = low - tm
    pct = round(diff / tm * 100, 3) if tm else ""
    st = "OK" if tm and abs(diff) / tm < 0.01 else ("HIGHER" if diff > 0 else "LOWER")
    if not tm: st = "NO_TABLE_VALUE"
    out.append(dict(sample=s, file=os.path.basename(p), seq_bp=seq, lowercase_bp=low,
                    upper_N=upN, lower_n=lowN, table_masked_bp=tm,
                    diff_bp=diff, diff_pct=pct, status=st))
    tally[st] = tally.get(st, 0) + 1
    print(f"  [{i}/{len(rows)}] {s:38s} {st}", flush=True)

cols = ["sample", "file", "seq_bp", "lowercase_bp", "upper_N", "lower_n",
        "table_masked_bp", "diff_bp", "diff_pct", "status"]
with open(A.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in out: w.writerow(r)

print("\n===== 汇总 =====")
for k2, v in sorted(tally.items()):
    print(f"  {k2:16s} {v}")
print(f"[写出] {A.out}  ({len(out)} 行)")
