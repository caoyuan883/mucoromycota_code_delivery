#!/usr/bin/env python3
import re, gzip
from pathlib import Path
from collections import defaultdict

import os, sys
# [release refactor] module root from the environment (historical: /home/py/Desktop/12_faa_12_gff_technical_validation)
if not os.environ.get("REF12_BASE"):
    sys.exit("ERROR: set REF12_BASE to the 12-genome technical-validation folder")
BASE = Path(os.environ["REF12_BASE"])
FAA_DIR = BASE / "12_faa"
GFF_DIR = BASE / "12_gff"
OUT_DIR = BASE / "12_faa_longest_isoform"
REPORT  = BASE / "longest_isoform_report.tsv"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def op(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "r")

def attrs(s):
    d = {}
    for kv in s.rstrip(";").split(";"):
        kv = kv.strip()
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def cands(val):
    """有序候选 key；丢弃门户名等歧义 token"""
    if not val: return []
    v = val.strip()
    out = [v, re.sub(r"\.\d+$", "", v)]
    for pre in ("cds-", "rna-", "gene-", "CDS:", "mRNA:"):
        if v.startswith(pre):
            t = v[len(pre):]
            out += [t, re.sub(r"\.\d+$", "", t)]
    if "|" in v:
        f = [x for x in v.split("|") if x]
        pick = ([f[-1]] if f else []) + [x for x in f if x.isdigit()]
        for c in pick:
            out += [c, re.sub(r"\.\d+$", "", c)]
    seen, res = set(), []
    for x in out:
        if x and x not in seen:
            seen.add(x); res.append(x)
    return res

def parse_gff(path):
    lines = []
    with op(path) as fh:
        for line in fh:
            if not line or line[0] == "#": continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9: lines.append((f[2], attrs(f[8])))
    mrna2gene, genes = {}, set()
    for t, a in lines:
        if t == "gene":
            if a.get("ID"): genes.add(a["ID"])
        elif t in ("mRNA", "transcript"):
            fid = a.get("ID", "")
            if fid: mrna2gene[fid] = a.get("Parent", "").split(",")[0] or fid
    pm, gm = defaultdict(set), defaultdict(set)
    for t, a in lines:
        if t in ("mRNA", "transcript"):
            fid = a.get("ID", "")
            gid = mrna2gene.get(fid, fid)
            pv = [a.get("proteinId",""), a.get("protein_id","")]
            gv = [fid, a.get("Name","")]
        elif t == "CDS":
            par = a.get("Parent", "").split(",")[0]
            gid = mrna2gene.get(par, par)
            pv = [a.get("protein_id",""), a.get("proteinId","")]
            gv = [a.get("ID",""), a.get("Name",""), par]
        else:
            continue
        if not gid: continue
        for v in pv:
            for k in cands(v): pm[k].add(gid)
        for v in gv:
            for k in cands(v): gm[k].add(gid)
    P = {k: next(iter(g)) for k, g in pm.items() if len(g) == 1}
    G = {k: next(iter(g)) for k, g in gm.items() if len(g) == 1}
    n_amb = sum(1 for g in pm.values() if len(g) > 1)
    return (P, G), len(genes), n_amb

def read_fasta(path):
    name, seq = None, []
    with op(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None: yield name, "".join(seq)
                name, seq = line[1:], []
            else: seq.append(line.strip())
    if name is not None: yield name, "".join(seq)

def find(d, stem, exts):
    for e in exts:
        p = d / (stem + e)
        if p.is_file(): return p
    return None

rows = []
for faa in sorted(p for p in FAA_DIR.iterdir()
                  if p.suffix in (".faa", ".fa", ".fasta")):
    stem = faa.stem
    gff = find(GFF_DIR, stem, [".gff3", ".gff", ".gff3.gz", ".gff.gz"])
    if gff is None:
        print(f"[缺GFF] {stem}"); continue
    key2gene, n_gene_gff, n_amb = parse_gff(gff)
    best, n_in, n_unmap, first_hdr = {}, 0, 0, None
    for hdr, seq in read_fasta(faa):
        n_in += 1
        if first_hdr is None: first_hdr = hdr.split()[0]
        s = seq.replace("*", "").replace("-", "").strip().upper()
        pid = hdr.split()[0]
        P, G = key2gene
        gid = None
        for k in cands(pid):
            if k in P: gid = P[k]; break
        if gid is None:
            for k in cands(pid):
                if k in G: gid = G[k]; break
        if gid is None:
            n_unmap += 1; gid = "__UNMAPPED__" + pid
        cur = best.get(gid)
        if cur is None or len(s) > cur[0] or (len(s) == cur[0] and pid < cur[2]):
            best[gid] = (len(s), hdr, pid, s)
    out = OUT_DIR / f"{stem}_longest_isoform.faa"
    with out.open("w") as fh:
        for gid in sorted(best, key=lambda g: best[g][2]):
            _, hdr, _, s = best[gid]
            fh.write(f">{hdr}\n")
            for i in range(0, len(s), 60): fh.write(s[i:i+60] + "\n")
    rate = 100.0 * (n_in - n_unmap) / n_in if n_in else 0.0
    flag = ""
    if rate < 99: flag += "  <-映射率低"
    if len(best) < 0.5 * n_gene_gff: flag += "  <-输出远少于GFF基因数"
    rows.append((stem, n_in, n_gene_gff, len(best), n_unmap,
                 n_amb, f"{rate:.2f}", first_hdr))
    print(f"{stem:32s} 输入{n_in:6d} GFF基因{n_gene_gff:6d} "
          f"输出{len(best):6d} 未映射{n_unmap:5d} 歧义key{n_amb:6d}{flag}")

with REPORT.open("w") as fh:
    fh.write("sample\tn_proteins_in\tn_genes_gff\tn_out\tn_unmapped\t"
             "n_ambiguous_keys\tmap_rate_pct\tfirst_fasta_id\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
print(f"\n共 {len(rows)} 株 -> {OUT_DIR}\n报告 -> {REPORT}")
