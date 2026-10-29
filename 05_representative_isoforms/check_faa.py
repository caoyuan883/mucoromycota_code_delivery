import sys, os, re, glob, hashlib, statistics
from collections import Counter
D = sys.argv[1]
OUT = os.path.dirname(os.path.abspath(D))
SUF = "_longest_isoform.faa"
STD = set("ACDEFGHIKLMNPQRSTVWY")
EXT = set("XBZJUO")

def norm(x):
    x = os.path.basename(x.strip().split()[0].split("\t")[0]) if x.strip() else ""
    for s in (SUF, ".faa", ".fa", ".fasta", "_longest_isoform"):
        if x.endswith(s): x = x[:-len(s)]
    return x

rows, warn = [], []
files = sorted(glob.glob(os.path.join(D, "*.faa")))
for f in files:
    b = os.path.basename(f)
    s = b[:-len(SUF)] if b.endswith(SUF) else b[:-4]
    if not b.endswith(SUF): warn.append(f"命名不符: {b}")
    if re.search(r"[^A-Za-z0-9_.\-]", s): warn.append(f"株名含特殊字符: {b}")
    md5 = hashlib.md5(open(f, "rb").read()).hexdigest()
    ids, lens, seq = [], [], []
    term = inter = nonstd = ext = empty = 0
    badline = 0
    def flush():
        global term, inter, nonstd, ext, empty
        q = "".join(seq).upper()
        if q.endswith("*"): term += 1; q = q[:-1]
        if "*" in q: inter += 1
        if not q: empty += 1
        c = set(q) - STD
        if c & EXT: ext += 1
        if c - EXT - {"*"}: nonstd += 1
        lens.append(len(q))
    with open(f) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if ids: flush()
                ids.append(line[1:].split()[0] if line[1:].strip() else ""); seq = []
            elif line:
                if not ids: badline += 1
                seq.append(line.strip())
        if ids: flush()
    dup = sum(v - 1 for v in Counter(ids).values() if v > 1)
    blank = sum(1 for i in ids if not i)
    n = len(ids)
    st = "OK"
    prob = []
    if n == 0: prob.append("无序列")
    if dup: prob.append(f"重复ID{dup}")
    if blank: prob.append(f"空ID{blank}")
    if empty: prob.append(f"空序列{empty}")
    if inter: prob.append(f"内部终止{inter}")
    if nonstd: prob.append(f"非法字符{nonstd}")
    if badline: prob.append("首行非header")
    if prob: st = ";".join(prob)
    rows.append([s, b, n, sum(lens),
                 min(lens) if lens else 0,
                 int(statistics.median(lens)) if lens else 0,
                 max(lens) if lens else 0,
                 term, inter, ext, dup, ids[0] if ids else "",
                 os.path.getsize(f), md5, st])

hdr = ["strain","file","n_proteins","total_aa","min_len","median_len","max_len",
       "terminal_stop","internal_stop","has_XBZJUO","dup_ids","first_id","bytes","md5","status"]
with open(os.path.join(OUT, "02_manifest.tsv"), "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")

strains = {r[0] for r in rows}
md5c = Counter(r[13] for r in rows)
rep = [f"目录: {D}", f"faa 文件数: {len(files)}",
       f"总蛋白数: {sum(r[2] for r in rows):,}",
       f"status 非 OK 的株: {sum(1 for r in rows if r[14] != 'OK')}",
       f"含末尾 * 的株: {sum(1 for r in rows if r[7] > 0)}"]
for m, c in md5c.items():
    if c > 1: rep.append("内容完全相同的文件: " + ", ".join(r[1] for r in rows if r[13] == m))
rep += warn
for r in rows:
    if r[14] != "OK": rep.append(f"问题: {r[1]}\t{r[14]}")
for t in sorted(glob.glob(os.path.join(D, "*.txt")) + glob.glob(os.path.join(OUT, "*.txt"))):
    lines = [l for l in open(t, errors="replace") if l.strip() and not l.startswith("#")]
    names = [norm(l) for l in lines]
    ns = set(names)
    rep.append(f"\n== 名单 {t}: 非空行 {len(lines)}，去重后 {len(ns)}")
    d = [k for k, v in Counter(names).items() if v > 1]
    if d: rep.append("  名单内重复: " + ", ".join(d))
    a = sorted(ns - strains); b = sorted(strains - ns)
    if a: rep.append("  名单有、文件没有: " + ", ".join(a))
    if b: rep.append("  文件有、名单没有: " + ", ".join(b))
    if not a and not b: rep.append("  与 faa 文件完全一致")
open(os.path.join(OUT, "02_check_report.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
