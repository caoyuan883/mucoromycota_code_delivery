#!/usr/bin/env python3
# omark_contam_v5.py —— 只读取, 不过滤、不修改任何原始文件
# 在 v4 (OMArk 标签 -> contig 分类) 基础上, 加入 contig 真实长度与 GC 作为独立证据
# 输出: ~/Desktop/omark_contam_v5
import os, re, csv, glob, gzip
from collections import defaultdict, Counter

HOME = os.path.expanduser("~")
SRC = os.path.join(HOME, "Desktop/omark_rerun_acaulospora")
GTF_DIR = os.path.join(HOME, "Desktop/Mucoromycota_BRAKER_longest_isoform/01_202_longest_isoform_gtf")
FA_DIR = os.path.join(HOME, "Desktop/10-fna")
FA_SUFFIX = "_genomic_softmasked.fa"
OUT = os.path.join(HOME, "Desktop/omark_contam_v5_acaulospora")
MIN_LEN_BASELINE = 5000     # 计算宿主 GC 基线时只用 >=5 kb 的 contig, 减少短片段噪声
DOMINANT = 0.8              # FA contig 上 污染/(污染+Consistent+Inconsistent) >= 0.8 视为"外源占绝对多数"
EXPECTED = ("mucor", "mortierell", "glomer", "fungi", "rhizop", "umbelops", "endogon",
            "diversispor", "gigaspor", "archaeospor", "paraglomer", "rhizophag")
for d in ("", "contig_level", "proposed_lists"):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)

def tsv(name, rows, fields=None):
    fields = fields or (list(rows[0].keys()) if rows else ["empty"])
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def opn(p): return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.endswith(".gz") else open(p, encoding="utf-8", errors="replace")

def parse_sum(p):
    clade, sp, on = "", [], False
    for s in opn(p).read().splitlines():
        m = re.search(r"selected clade was\s*(.+)$", s, re.I)
        if m: clade = m.group(1).strip()
        if "detected species" in s.lower(): on = True; continue
        if on:
            if s.startswith("#"): continue
            if not s.strip():
                if sp: on = False
                continue
            q = s.split("\t")
            if len(q) >= 4: sp.append(q[0].strip())
    return clade, sp

def parse_ump(p):
    cat, cur = {}, None
    for s in opn(p):
        s = s.strip()
        if not s: continue
        if s.startswith(">"):
            h = s[1:].strip().lower()
            cur = ("contam" if h.startswith("contamination") else "consistent" if h.startswith("consistent") else
                   "inconsistent" if h.startswith("inconsistent") else "unknown" if h.startswith("unknown") else None)
        elif cur:
            cat[s.split()[0]] = cur
    return cat

def gtf_maps(p):
    t2g, g2c = {}, {}
    for line in opn(p):
        if line.startswith("#"): continue
        q = line.rstrip("\n").split("\t")
        if len(q) < 9: continue
        c, feat, a = q[0], q[2], q[8].strip()
        mt = re.search(r'transcript_id "([^"]+)"', a); mg = re.search(r'gene_id "([^"]+)"', a)
        if mt:
            t = mt.group(1); g = mg.group(1) if mg else t.rsplit(".", 1)[0]
            t2g[t] = g; g2c[g] = c
        elif a and " " not in a and ";" not in a:
            if feat in ("transcript", "mRNA"): t2g[a] = a.rsplit(".", 1)[0]; g2c[t2g[a]] = c
            elif feat == "gene": g2c[a] = c
    return t2g, g2c

def to_gene(pid, t2g, g2c):
    cands = [pid, pid.split("|")[-1]]
    m = re.search(r"(g\d+\.t\d+)", pid)
    if m: cands.append(m.group(1))
    for x in cands:
        if x in t2g: return t2g[x]
        if x in g2c: return x
        if x.rsplit(".", 1)[0] in g2c: return x.rsplit(".", 1)[0]
    return None

def fasta_stats(p):
    st, name, L, GC, N = {}, None, 0, 0, 0
    def flush():
        if name is not None:
            acgt = L - N
            st[name] = (L, round(100 * GC / acgt, 2) if acgt else None)
    for line in opn(p):
        if line.startswith(">"):
            flush(); name = line[1:].split()[0]; L = GC = N = 0
        else:
            s = line.strip().upper()
            L += len(s); GC += s.count("G") + s.count("C")
            N += len(s) - (s.count("A") + s.count("C") + s.count("G") + s.count("T"))
    flush()
    return st

def wq(vals, q):
    """长度加权分位数; vals=[(gc, len)]"""
    v = sorted(x for x in vals if x[0] is not None)
    if not v: return None
    tot, acc = sum(l for _, l in v), 0
    for g, l in v:
        acc += l
        if acc >= q * tot: return g
    return v[-1][0]

def norm(x): return re.sub(r"[^a-z0-9]", "", x.lower())

ump = {os.path.basename(p)[:-4]: p for p in glob.glob(os.path.join(SRC, "**", "*.ump"), recursive=True)}
sm = {os.path.basename(p)[:-4]: p for p in glob.glob(os.path.join(SRC, "**", "*.sum"), recursive=True)}
gidx = {norm(re.sub(r"(_longest_isoform)?\.gtf(\.gz)?$", "", os.path.basename(p))): p
        for p in glob.glob(os.path.join(GTF_DIR, "**", "*.gtf*"), recursive=True)}
fidx = {norm(os.path.basename(p).split(FA_SUFFIX)[0]): p
        for p in glob.glob(os.path.join(FA_DIR, "*" + FA_SUFFIX + "*"))}
print(f"FASTA 文件: {len(fidx)}  (目录 {FA_DIR})")

summary, hist = [], []
for s in sorted(set(ump) & set(sm)):
    cat = parse_ump(ump[s])
    if not any(v == "contam" for v in cat.values()): continue
    gtf, fa = gidx.get(norm(s)), fidx.get(norm(s))
    if not gtf or not fa:
        print(f"[WARN] {s}: GTF={'OK' if gtf else '缺'} FASTA={'OK' if fa else '缺'}，跳过"); continue
    clade, sp = parse_sum(sm[s])
    suspect = not any(e in clade.lower() for e in EXPECTED)
    t2g, g2c = gtf_maps(gtf)
    st = fasta_stats(fa)
    miss_ctg = len({c for c in g2c.values() if c not in st})
    gcat = {}
    for pid, c in cat.items():
        g = to_gene(pid, t2g, g2c)
        if g: gcat[g] = c
    per = defaultdict(Counter)
    for g, c in g2c.items(): per[c][gcat.get(g, "not_in_omark")] += 1

    # 宿主基线: 含 Consistent 且无污染基因的 contig (>= MIN_LEN_BASELINE)
    host = [(st[c][1], st[c][0]) for c, k in per.items()
            if c in st and k["consistent"] > 0 and k["contam"] == 0 and st[c][0] >= MIN_LEN_BASELINE]
    if len(host) < 10:
        host = [(st[c][1], st[c][0]) for c, k in per.items() if c in st and k["consistent"] > 0 and k["contam"] == 0]
    h50, h05, h95 = wq(host, .5), wq(host, .05), wq(host, .95)

    rows = []
    for c, k in per.items():
        if k["contam"] == 0: continue
        L, gc = st.get(c, (None, None))
        cls = "FA" if k["consistent"] > 0 else ("FI" if k["inconsistent"] > 0 else "FO")
        assigned = k["contam"] + k["consistent"] + k["inconsistent"]
        share = k["contam"] / assigned if assigned else 0
        rows.append(dict(contig=c, length=L, gc=gc, n_genes=sum(k.values()), consistent=k["consistent"],
                         inconsistent=k["inconsistent"], contam=k["contam"], unknown=k["unknown"],
                         foreign_share=round(share, 3), cls=cls))
    fo = [(r["gc"], r["length"]) for r in rows if r["cls"] == "FO" and r["length"]]
    f50 = wq(fo, .5)

    def closer_foreign(gc):
        if gc is None or h50 is None or f50 is None: return None
        return abs(gc - f50) < abs(gc - h50)
    def outside_host(gc):
        if gc is None or h05 is None: return None
        return gc < h05 or gc > h95

    for r in rows:
        cf, oh = closer_foreign(r["gc"]), outside_host(r["gc"])
        r["gc_closer_to_foreign"], r["gc_outside_host_P5_P95"] = cf, oh
        r["fa_dominant"] = r["cls"] == "FA" and r["foreign_share"] >= DOMINANT
        if r["cls"] == "FO":
            r["proposal"] = "remove(FO)"
        elif r["fa_dominant"] and cf and oh:
            r["proposal"] = "remove_candidate(FA_dominant+GC)"
        elif r["fa_dominant"]:
            r["proposal"] = "keep_flag(FA_dominant, GC 不支持)"
        else:
            r["proposal"] = "keep_flag"
        if suspect: r["proposal"] = "HOLD(主体谱系可疑)"
    rows.sort(key=lambda r: (r["proposal"], -r["contam"]))
    tsv(f"contig_level/{s}.tsv", rows)
    for tag, cond in (("remove_FO", lambda r: r["proposal"] == "remove(FO)"),
                      ("remove_candidate_FA_GC", lambda r: r["proposal"].startswith("remove_candidate"))):
        with open(os.path.join(OUT, "proposed_lists", f"{s}.{tag}.contigs.txt"), "w") as f:
            f.write("".join(r["contig"] + "\n" for r in rows if cond(r)))

    # GC 直方图 (按 bp), 便于后续画图
    lab = {r["contig"]: r["cls"] for r in rows}
    for c, (L, gc) in st.items():
        if gc is None: continue
        k = per.get(c)
        grp = lab.get(c) or ("host_clean" if k and k["consistent"] > 0 else ("no_gene" if not k else "other"))
        hist.append((s, grp, int(gc // 1), L))

    tot_genes, tot_c = len(g2c), sum(1 for v in gcat.values() if v == "contam")
    asm_bp = sum(v[0] for v in st.values())
    def S(cond, key): return sum((r[key] or 0) for r in rows if cond(r))
    isFO = lambda r: r["cls"] == "FO"
    isUP = lambda r: r["proposal"].startswith("remove_candidate")
    isFAD = lambda r: r["fa_dominant"]
    fo_rows = [r for r in rows if isFO(r)]
    fo_cf = [r for r in fo_rows if r["gc_closer_to_foreign"]]
    summary.append(dict(
        sample=s, flag="HOLD(主体谱系可疑)" if suspect else "", selected_clade=clade, main=sp[0] if sp else "",
        assembly_bp=asm_bp, n_contigs=len(st), gtf_contigs_missing_in_fasta=miss_ctg,
        host_gc_median=h50, host_gc_P5=h05, host_gc_P95=h95, n_host_baseline_contigs=len(host),
        FO_gc_median=f50, gc_gap=round(f50 - h50, 2) if (f50 is not None and h50 is not None) else None,
        FO_contigs=len(fo_rows), FO_bp=S(isFO, "length"),
        FO_pct_contigs_gc_closer_foreign=round(100 * len(fo_cf) / len(fo_rows), 1) if fo_rows else None,
        FA_dominant_contigs=sum(1 for r in rows if isFAD(r)), FA_dominant_bp=S(isFAD, "length"),
        FA_dominant_contam=S(isFAD, "contam"),
        FA_upgrade_contigs=sum(1 for r in rows if isUP(r)), FA_upgrade_bp=S(isUP, "length"),
        FA_upgrade_contam=S(isUP, "contam"), FA_upgrade_consistent_lost=S(isUP, "consistent"),
        contam_total=tot_c,
        contam_removed_FO=S(isFO, "contam"), contam_removed_FO_plus_upgrade=S(isFO, "contam") + S(isUP, "contam"),
        pct_contam_removed_FO_plus_upgrade=round(100 * (S(isFO, "contam") + S(isUP, "contam")) / tot_c, 1) if tot_c else 0,
        genes_removed_FO_plus_upgrade=S(isFO, "n_genes") + S(isUP, "n_genes"),
        pct_genes_removed=round(100 * (S(isFO, "n_genes") + S(isUP, "n_genes")) / tot_genes, 2) if tot_genes else 0,
        bp_removed_FO_plus_upgrade=S(isFO, "length") + S(isUP, "length"),
        pct_bp_removed=round(100 * (S(isFO, "length") + S(isUP, "length")) / asm_bp, 2) if asm_bp else 0,
    ))

tsv("01_sample_summary.tsv", summary)
agg = Counter()
for s, g, b, L in hist: agg[(s, g, b)] += L
tsv("02_gc_histogram_bp.tsv", [dict(sample=s, group=g, gc_bin=b, bp=v) for (s, g, b), v in sorted(agg.items())])

print("\nsample\thostGC(P5-P95)\tFO_GC\tFO_closer%\tFAdom(n/contam)\tupgrade(n/contam/consLost)\tremoved contam%\tgenes%\tbp%\tflag")
for r in summary:
    print(f"{r['sample']}\t{r['host_gc_median']}({r['host_gc_P5']}-{r['host_gc_P95']})\t{r['FO_gc_median']}\t"
          f"{r['FO_pct_contigs_gc_closer_foreign']}\t{r['FA_dominant_contigs']}/{r['FA_dominant_contam']}\t"
          f"{r['FA_upgrade_contigs']}/{r['FA_upgrade_contam']}/{r['FA_upgrade_consistent_lost']}\t"
          f"{r['pct_contam_removed_FO_plus_upgrade']}\t{r['pct_genes_removed']}\t{r['pct_bp_removed']}\t{r['flag']}")
    if r["gtf_contigs_missing_in_fasta"]:
        print(f"  [WARN] {r['sample']}: {r['gtf_contigs_missing_in_fasta']} 个 GTF contig 名在 FASTA 中找不到")
print(f"\n未过滤任何文件。输出: {OUT}")
