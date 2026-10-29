#!/usr/bin/env python3
# omark_filter_v6.py —— 按最终规则生成"过滤版"文件; 原始文件只读, 一律不改
# 规则 RU: 移除 contig 当且仅当 [FO 或 (FA 且 外源占比>=0.8 且 GC 更近外源)] 且 contig GC 在宿主 P5–P95 之外
# 例外 (仅 Acaulospora_colombiana): 另移除装配中所有 GC >= 43% 的 contig (双峰中的担子菌组分)
# 输出: ~/Desktop/omark_filtered_v6
import os, re, csv, glob, gzip, hashlib
from collections import defaultdict

HOME = os.path.expanduser("~")
V5_DIRS = [os.path.join(HOME, "Desktop/omark_contam_v5_acaulospora/contig_level"),   # 优先: Acaulospora 新结果
           os.path.join(HOME, "Desktop/omark_contam_v5/contig_level")]
GTF_DIR = os.path.join(HOME, "Desktop/Mucoromycota_BRAKER_longest_isoform/01_202_longest_isoform_gtf")
FAA_DIRS = [os.path.join(HOME, "Desktop/Mucoromycota_BRAKER_longest_isoform/02_longest_isoform_proteins"),
            os.path.join(HOME, "Desktop/Mucoromycota_BRAKER_longest_isoform")]
FA_DIR = os.path.join(HOME, "Desktop/10-fna")
FA_SUFFIX = "_genomic_softmasked.fa"
OUT = os.path.join(HOME, "Desktop/omark_filtered_v6")
GC_EXCEPTION = {"Acaulospora_colombiana": 43.0}
DOMINANT = 0.8

def opn(p): return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.endswith(".gz") else open(p, encoding="utf-8", errors="replace")
def truthy(x): return str(x).strip().lower() == "true"
def tsv(path, rows, fields=None):
    fields = fields or (list(rows[0].keys()) if rows else ["empty"])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

# ---- 1. 读取 v5 contig 表 (Acaulospora 用新版) ----
v5 = {}
for d in V5_DIRS:
    for p in glob.glob(os.path.join(d, "*.tsv")):
        s = os.path.basename(p)[:-4]
        v5.setdefault(s, p)
if not v5: raise SystemExit("[ERROR] 找不到 v5 contig_level 表")
print(f"v5 contig 表: {len(v5)} 株  (Acaulospora 来源: {v5.get('Acaulospora_colombiana','缺')})")

def find(dirs, sample, pats):
    for d in dirs:
        for pat in pats:
            hits = glob.glob(os.path.join(d, "**", pat.format(s=sample)), recursive=True)
            if hits: return sorted(hits)[0]
    return None

summary, manifest = [], []
for s in sorted(v5):
    rows = list(csv.DictReader(open(v5[s], encoding="utf-8"), delimiter="\t"))
    gtf = find([GTF_DIR], s, ["{s}_longest_isoform.gtf", "{s}*.gtf"])
    fa = find([FA_DIR], s, ["{s}" + FA_SUFFIX, "{s}*.fa", "{s}*.fna", "{s}*.fasta"])
    faa = find(FAA_DIRS, s, ["{s}_longest_isoform.faa", "{s}*longest*.faa"])
    if not gtf or not fa:
        print(f"[WARN] {s}: GTF={'OK' if gtf else '缺'} FASTA={'OK' if fa else '缺'}, 跳过"); continue
    od = os.path.join(OUT, s); os.makedirs(od, exist_ok=True)

    # ---- 2. 规则 RU ----
    reason = {}
    for r in rows:
        fo = r["cls"] == "FO"
        fad = r["cls"] == "FA" and float(r["foreign_share"]) >= DOMINANT and truthy(r["gc_closer_to_foreign"])
        if (fo or fad) and truthy(r["gc_outside_host_P5_P95"]):
            reason[r["contig"]] = "RU_FO" if fo else "RU_FA_dominant"
    contam_on = {r["contig"]: int(r["contam"]) for r in rows}

    # ---- 3. 读装配, 计算 GC, 应用例外 ----
    seqs_len, seqs_gc, order = {}, {}, []
    name, L, GC, N = None, 0, 0, 0
    def flush():
        if name is not None:
            seqs_len[name] = L; acgt = L - N
            seqs_gc[name] = 100 * GC / acgt if acgt else None; order.append(name)
    for line in opn(fa):
        if line.startswith(">"):
            flush(); name = line[1:].split()[0]; L = GC = N = 0
        else:
            t = line.strip().upper(); L += len(t); GC += t.count("G") + t.count("C")
            N += len(t) - sum(t.count(b) for b in "ACGT")
    flush()
    thr = GC_EXCEPTION.get(s)
    if thr is not None:
        for c in order:
            if seqs_gc[c] is not None and seqs_gc[c] >= thr and c not in reason:
                reason[c] = f"EXC_GC>={thr:g}"

    # ---- 4. GTF: gene/transcript -> contig; 写过滤版 GTF ----
    t2g, g2c = {}, {}
    gtf_out = os.path.join(od, f"{s}_longest_isoform.filtered.gtf")
    kept_lines = dropped_lines = 0
    with opn(gtf) as fi, open(gtf_out, "w", encoding="utf-8") as fo_:
        for line in fi:
            if line.startswith("#"): fo_.write(line); continue
            q = line.rstrip("\n").split("\t")
            if len(q) < 9: fo_.write(line); continue
            c, feat, a = q[0], q[2], q[8].strip()
            mt = re.search(r'transcript_id "([^"]+)"', a); mg = re.search(r'gene_id "([^"]+)"', a)
            if mt:
                t = mt.group(1); g = mg.group(1) if mg else t.rsplit(".", 1)[0]; t2g[t] = g; g2c[g] = c
            elif a and " " not in a and ";" not in a:
                if feat in ("transcript", "mRNA"): t2g[a] = a.rsplit(".", 1)[0]; g2c[t2g[a]] = c
                elif feat == "gene": g2c[a] = c
            if c in reason: dropped_lines += 1
            else: fo_.write(line); kept_lines += 1
    rm_genes = sorted(g for g, c in g2c.items() if c in reason)
    with open(os.path.join(od, f"{s}.removed_genes.txt"), "w") as f:
        f.write("".join(g + "\n" for g in rm_genes))

    # ---- 5. 过滤版 FASTA + 移除序列 FASTA ----
    fa_keep = os.path.join(od, f"{s}_genomic_softmasked.filtered.fa")
    fa_rm = os.path.join(od, f"{s}.removed_contigs.fa")
    with opn(fa) as fi, open(fa_keep, "w") as fk, open(fa_rm, "w") as fr:
        out = fk
        for line in fi:
            if line.startswith(">"):
                out = fr if line[1:].split()[0] in reason else fk
            out.write(line)

    # ---- 6. 过滤版蛋白 FAA (如有) ----
    n_prot_in = n_prot_out = None
    faa_out = None
    if faa:
        faa_out = os.path.join(od, f"{s}_longest_isoform.filtered.faa")
        n_prot_in = n_prot_out = 0; keep = True
        with opn(faa) as fi, open(faa_out, "w") as fo_:
            for line in fi:
                if line.startswith(">"):
                    pid = line[1:].split()[0]; n_prot_in += 1
                    g = t2g.get(pid) or (pid if pid in g2c else pid.rsplit(".", 1)[0])
                    keep = g2c.get(g) not in reason
                    n_prot_out += keep
                if keep: fo_.write(line)

    # ---- 7. 审计表 ----
    audit = []
    for c in sorted(reason, key=lambda x: (reason[x], x)):
        audit.append(dict(contig=c, reason=reason[c], length=seqs_len.get(c), gc=round(seqs_gc[c], 2) if seqs_gc.get(c) is not None else None,
                          contam_genes=contam_on.get(c, 0), genes=sum(1 for g in rm_genes if g2c[g] == c)))
    tsv(os.path.join(od, f"{s}.removed_contigs_audit.tsv"), audit)
    asm_bp = sum(seqs_len.values()); rm_bp = sum(seqs_len.get(c, 0) for c in reason)
    contam_total = sum(contam_on.values()); contam_rm = sum(contam_on.get(c, 0) for c in reason)
    by = defaultdict(int)
    for c in reason: by[reason[c]] += 1
    summary.append(dict(
        sample=s, rule="RU" + (f"+GC>={thr:g}" if thr else ""),
        contigs_before=len(order), contigs_removed=len(reason),
        removed_RU_FO=by["RU_FO"], removed_RU_FA_dominant=by["RU_FA_dominant"],
        removed_GC_exception=sum(v for k, v in by.items() if k.startswith("EXC")),
        bp_before=asm_bp, bp_removed=rm_bp, pct_bp_removed=round(100 * rm_bp / asm_bp, 2),
        genes_before=len(g2c), genes_removed=len(rm_genes), pct_genes_removed=round(100 * len(rm_genes) / len(g2c), 2),
        contam_genes_total=contam_total, contam_genes_removed=contam_rm,
        pct_contam_removed=round(100 * contam_rm / contam_total, 1) if contam_total else 0,
        proteins_in=n_prot_in, proteins_out=n_prot_out,
        gtf_lines_kept=kept_lines, gtf_lines_dropped=dropped_lines,
        faa_status="OK" if faa else "未找到原蛋白 FAA (只生成 GTF/FASTA/基因清单)"))
    for p in [gtf_out, fa_keep, fa_rm, faa_out, os.path.join(od, f"{s}.removed_genes.txt"), os.path.join(od, f"{s}.removed_contigs_audit.tsv")]:
        if p and os.path.exists(p):
            manifest.append(dict(sample=s, file=os.path.relpath(p, OUT), bytes=os.path.getsize(p), sha256=sha256(p)))
    for src in (gtf, fa, faa):
        if src: manifest.append(dict(sample=s, file="SOURCE:" + src, bytes=os.path.getsize(src), sha256=sha256(src)))

tsv(os.path.join(OUT, "00_filter_summary.tsv"), summary)
tsv(os.path.join(OUT, "MANIFEST.tsv"), manifest, ["sample", "file", "bytes", "sha256"])
print("\nsample\trule\tcontigs_rm(FO/FAd/GCexc)\tbp_rm%\tgenes_rm(%)\tcontam_rm%\tproteins in->out")
for r in summary:
    print(f"{r['sample']}\t{r['rule']}\t{r['contigs_removed']}({r['removed_RU_FO']}/{r['removed_RU_FA_dominant']}/{r['removed_GC_exception']})\t"
          f"{r['pct_bp_removed']}\t{r['genes_removed']}({r['pct_genes_removed']}%)\t{r['pct_contam_removed']}\t{r['proteins_in']}->{r['proteins_out']}")
    if r["faa_status"] != "OK": print(f"  [NOTE] {r['faa_status']}")
print(f"\n原始文件未改动。输出: {OUT}")
