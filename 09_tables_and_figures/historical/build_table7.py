#!/usr/bin/env python3
"""表7 提取: 从 longest_isoform GTF/FAA 统计基因结构指标。
用法:
  build_table7.py --gtf-dir DIR --faa-dir DIR [--raw-gtf-dir DIR] [--out FILE]
"""
import os, re, sys, csv, argparse
from statistics import median
from collections import defaultdict

RE_TID = re.compile(r'transcript_id\s+"([^"]+)"')
RE_GID = re.compile(r'gene_id\s+"([^"]+)"')

def ids(feat, attr):
    t = RE_TID.search(attr); g = RE_GID.search(attr)
    tid = t.group(1) if t else None
    gid = g.group(1) if g else None
    if tid is None and feat == "transcript":
        tid = attr.strip().strip(';').strip('"')
    if gid is None and feat == "gene":
        gid = attr.strip().strip(';').strip('"')
    if gid is None and tid:
        gid = re.sub(r'\.t\d+$', '', tid)
    return gid, tid

def parse_gtf(path):
    genes, tx = set(), set()
    exon = defaultdict(int); cds = defaultdict(int)
    feats = defaultdict(int)
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            feat, attr = f[2], f[8]
            feats[feat] += 1
            gid, tid = ids(feat, attr)
            if gid: genes.add(gid)
            if tid: tx.add(tid)
            if feat == "exon" and tid: exon[tid] += 1
            elif feat == "CDS" and tid: cds[tid] += 1
    return genes, tx, exon, cds, feats

def count_faa(path):
    n = 0; total = 0
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"): n += 1
            else: total += len(line.strip().replace("*", ""))
    return n, total

def stem(fn, suf):
    b = os.path.basename(fn)
    for s in suf:
        if b.endswith(s): return b[:-len(s)]
    return os.path.splitext(b)[0]

a = argparse.ArgumentParser()
a.add_argument("--gtf-dir", required=True)
a.add_argument("--faa-dir", required=True)
a.add_argument("--raw-gtf-dir")
a.add_argument("--out", default="table7_gene_structure.tsv")
A = a.parse_args()
P = os.path.expanduser

gdir, fdir = P(A.gtf_dir), P(A.faa_dir)
gtfs = {stem(f, ["_longest_isoform.gtf", ".gtf"]): os.path.join(gdir, f)
        for f in sorted(os.listdir(gdir)) if f.endswith(".gtf")}
faas = {stem(f, ["_longest_isoform.faa", ".faa", ".fa", ".fasta"]): os.path.join(fdir, f)
        for f in sorted(os.listdir(fdir)) if f.endswith((".faa", ".fa", ".fasta"))}
raws = {}
if A.raw_gtf_dir:
    rd = P(A.raw_gtf_dir)
    raws = {stem(f, ["_braker.gtf", ".gtf"]): os.path.join(rd, f)
            for f in sorted(os.listdir(rd)) if f.endswith(".gtf")}
print(f"GTF {len(gtfs)} | FAA {len(faas)} | RAW {len(raws)}")

rows = []
for i, (s, gp) in enumerate(sorted(gtfs.items()), 1):
    genes, tx, exon, cds, feats = parse_gtf(gp)
    use = "exon" if exon and len(exon) >= len(tx) else "CDS"
    per = exon if use == "exon" else cds
    counts = [per.get(t, 0) for t in tx if per.get(t, 0) > 0]
    mono = sum(1 for c in counts if c == 1)
    mult = sum(1 for c in counts if c > 1)
    mm = [c for c in counts if c > 1]
    fp = faas.get(s) or next((v for k, v in faas.items()
                              if k.startswith(s) or s.startswith(k)), None)
    nfaa, aatot = count_faa(fp) if fp else (None, None)
    rg = raws.get(s) or next((v for k, v in raws.items()
                              if k.startswith(s) or s.startswith(k)), None)
    if rg:
        rgen, rtx, _, _, _ = parse_gtf(rg)
        n_rg, n_rt = len(rgen), len(rtx)
    else:
        n_rg = n_rt = ""
    rows.append(dict(
        Species=s,
        n_genes=len(genes), n_transcripts=len(tx),
        n_proteins_faa=nfaa if nfaa is not None else "",
        exon_feature_used=use,
        n_mono_exon=mono, n_multi_exon=mult,
        mono_mult_ratio=round(mono / mult, 3) if mult else "",
        median_exons_multiexon=int(median(mm)) if mm else "",
        max_exons=max(counts) if counts else "",
        mean_protein_len=round(aatot / nfaa, 1) if nfaa else "",
        n_genes_raw=n_rg, n_transcripts_raw=n_rt,
        isoforms_per_gene=round(n_rt / n_rg, 3) if rg and n_rg else "",
        check_genes_eq_faa="OK" if nfaa == len(genes) else "MISMATCH",
        check_tx_eq_genes="OK" if len(tx) == len(genes) else "MISMATCH",
        gtf_features=";".join(f"{k}:{v}" for k, v in sorted(feats.items())),
    ))
    print(f"  [{i}/{len(gtfs)}] {s:40s} genes={len(genes)} tx={len(tx)} {use}", flush=True)

cols = list(rows[0].keys())
with open(A.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows: w.writerow(r)

nb1 = sum(1 for r in rows if r["check_genes_eq_faa"] != "OK")
nb2 = sum(1 for r in rows if r["check_tx_eq_genes"] != "OK")
print(f"\n===== 汇总 =====")
print(f"  样本数 {len(rows)}")
print(f"  基因数 != FAA 条数 : {nb1}")
print(f"  转录本数 != 基因数 : {nb2}")
print(f"  外显子特征使用情况 : "
      + str({u: sum(1 for r in rows if r['exon_feature_used'] == u) for u in {r['exon_feature_used'] for r in rows}}))
print(f"[写出] {A.out} ({len(rows)} 行 x {len(cols)} 列)")
