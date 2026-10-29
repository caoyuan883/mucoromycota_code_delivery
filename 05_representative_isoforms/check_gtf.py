import sys, os, re, glob, hashlib, statistics
from collections import Counter, defaultdict
G, P = sys.argv[1], sys.argv[2]
OUT = os.path.dirname(os.path.abspath(G))
SUF = "_longest_isoform"
rt = re.compile(r'transcript_id "([^"]+)"'); rg = re.compile(r'gene_id "([^"]+)"')

def strain(f):
    s = os.path.basename(f).rsplit(".", 1)[0]
    return s[:-len(SUF)] if s.endswith(SUF) else s

def faa_ids(s):
    f = os.path.join(P, s + SUF + ".faa")
    if not os.path.exists(f): return None
    return {l[1:].split()[0] for l in open(f) if l.startswith(">")}

rows, rep = [], []
files = sorted(glob.glob(os.path.join(G, "*.gtf")))
for f in files:
    s = strain(f); b = os.path.basename(f)
    md5 = hashlib.md5(open(f, "rb").read()).hexdigest()
    feat = Counter(); seqids = set(); genes = set(); g2t = defaultdict(set)
    cds = Counter(); exon = Counter(); sc = set(); ec = set(); txl = set()
    bad = coord = strand = noid = 0
    for line in open(f):
        if line.startswith("#") or not line.strip(): continue
        c = line.rstrip("\n").split("\t")
        if len(c) != 9: bad += 1; continue
        seqids.add(c[0]); ft = c[2]; feat[ft] += 1
        try:
            if int(c[3]) > int(c[4]) or int(c[3]) < 1: coord += 1
        except ValueError: coord += 1
        if c[6] not in "+-": strand += 1
        a = c[8]; t = rt.search(a); g = rg.search(a)
        tid = t.group(1) if t else None; gid = g.group(1) if g else None
        if ft == "gene": genes.add(gid or a.strip()); continue
        if ft in ("transcript", "mRNA") and not tid: tid = a.strip()
        if not tid: noid += 1; continue
        if gid: g2t[gid].add(tid); genes.add(gid)
        if ft in ("transcript", "mRNA"): txl.add(tid)
        elif ft == "CDS": cds[tid] += 1
        elif ft == "exon": exon[tid] += 1
        elif ft == "start_codon": sc.add(tid)
        elif ft == "stop_codon": ec.add(tid)
    tx = set(cds) | txl | set(exon)
    multi_iso = sum(1 for v in g2t.values() if len(v) > 1)
    nocds = len(tx - set(cds))
    mono = sum(1 for t in cds if cds[t] == 1); multi = len(cds) - mono
    ratio = round(mono / multi, 3) if multi else "NA"
    nostart = len(set(cds) - sc); nostop = len(set(cds) - ec)
    fi = faa_ids(s)
    if fi is None: fo = go = "NA"; fstat = "无对应faa"
    else:
        fo = len(fi - tx); go = len(tx - fi); fstat = "一致" if fo == 0 and go == 0 else "不一致"
    prob = []
    if not tx: prob.append("无转录本")
    if bad: prob.append(f"非9列{bad}")
    if coord: prob.append(f"坐标异常{coord}")
    if strand: prob.append(f"链异常{strand}")
    if noid: prob.append(f"缺transcript_id{noid}")
    if multi_iso: prob.append(f"多isoform基因{multi_iso}")
    if nocds: prob.append(f"无CDS转录本{nocds}")
    if fstat != "一致": prob.append(f"faa{fstat}")
    rows.append([s, b, len(genes), len(tx), len(cds), feat["CDS"], mono, multi, ratio,
                 nostart, nostop, len(seqids), fo, go, os.path.getsize(f), md5,
                 ";".join(prob) if prob else "OK"])

hdr = ["strain","file","n_genes","n_transcripts","n_tx_with_CDS","n_CDS_lines","mono_exon","multi_exon",
       "mono_multi_ratio","tx_no_start_codon","tx_no_stop_codon","n_seqids",
       "faa_only_ids","gtf_only_ids","bytes","md5","status"]
with open(os.path.join(OUT, "01_gtf_manifest.txt"), "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")

gs = {r[0] for r in rows}
ps = {os.path.basename(f)[:-len(SUF + ".faa")] for f in glob.glob(os.path.join(P, "*" + SUF + ".faa"))}
rat = [r[8] for r in rows if r[8] != "NA"]
rep += [f"GTF 目录: {G}", f"GTF 文件数: {len(files)}", f"faa 文件数: {len(ps)}",
        f"总基因数: {sum(r[2] for r in rows):,}", f"总转录本数: {sum(r[3] for r in rows):,}",
        f"单/多外显子比 中位数: {statistics.median(rat) if rat else 'NA'}",
        f"status 非 OK 的株: {sum(1 for r in rows if r[16] != 'OK')}",
        f"含多 isoform 基因的株: {sum(1 for r in rows if 'isoform' in r[16])}",
        f"与 faa ID 完全一致的株: {sum(1 for r in rows if r[12] == 0 and r[13] == 0)}"]
a = sorted(gs - ps); b = sorted(ps - gs)
if a: rep.append("有 GTF 无 faa: " + ", ".join(a))
if b: rep.append("有 faa 无 GTF: " + ", ".join(b))
for m, c in Counter(r[15] for r in rows).items():
    if c > 1: rep.append("内容完全相同的 GTF: " + ", ".join(r[1] for r in rows if r[15] == m))
for r in rows:
    if r[16] != "OK": rep.append(f"问题: {r[1]}\t{r[16]}")
open(os.path.join(OUT, "01_gtf_check_report.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
