#!/usr/bin/env python3
import os, re, sys, argparse
from collections import OrderedDict, defaultdict

TBL_ROW = re.compile(r"^\s*([A-Za-z][^:]*?):?\s+(\d+)\s+(\d+)\s*bp\s+([\d.]+)\s*%\s*$")
KEEP = OrderedDict([("Retroelements","retroelement"),("SINEs","SINE"),("LINEs","LINE"),
    ("LTR elements","LTR"),("DNA transposons","DNA_TE"),("Rolling-circles","RC_Helitron"),
    ("Unclassified","unclassified"),("Small RNA","small_RNA"),("Satellites","satellite"),
    ("Simple repeats","simple_repeat"),("Low complexity","low_complexity")])

def rd(p):
    return open(p, "r", errors="replace").read()

def g1(pat, txt, flags=re.M):
    m = re.search(pat, txt, flags)
    return m.group(1).strip() if m else None

def norm(fn):
    b = os.path.basename(fn)
    for s in (".fna.tbl",".fa.tbl",".tbl","-rmod.log","_rmod.log","-families.fa",
              "_families.fa",".fna.out",".out",".log",".txt"):
        if b.endswith(s):
            b = b[:-len(s)]; break
    return re.sub(r"[._-]+$", "", re.sub(r"_genomic$", "", b))

def index(d, exts):
    r = {}
    if d and os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            if any(fn.endswith(e) for e in exts):
                r.setdefault(norm(fn), os.path.join(d, fn))
    return r

def look(idx, s):
    if s in idx: return idx[s]
    for k, v in idx.items():
        if k.startswith(s) or s.startswith(k): return v
    f = re.sub(r"[^a-z0-9]", "", s.lower())
    for k, v in idx.items():
        if re.sub(r"[^a-z0-9]", "", k.lower()) == f: return v
    return None

def p_tbl(p):
    t = rd(p); d = {}
    d["sequences"] = g1(r"sequences:\s+(\d+)", t)
    d["total_length_bp"] = g1(r"total length:\s+(\d+)\s*bp", t)
    d["GC_percent"] = g1(r"GC level:\s+([\d.]+)\s*%", t)
    m = re.search(r"bases masked:\s+(\d+)\s*bp\s*\(\s*([\d.]+)\s*%\)", t)
    d["total_masked_bp"], d["total_masked_pct"] = (m.group(1), m.group(2)) if m else (None, None)
    m = re.search(r"Total interspersed repeats:\s+(\d+)\s*bp\s+([\d.]+)\s*%", t)
    d["interspersed_bp"], d["interspersed_pct"] = (m.group(1), m.group(2)) if m else (None, None)
    d["rm_library"] = (g1(r"RepeatMasker Combined Database:\s*(\S+)", t)
                       or g1(r'compared to classified sequences in\s*"?([^"\n]+)"?', t))
    rows = {}
    for ln in t.splitlines():
        m = TBL_ROW.match(ln)
        if m: rows.setdefault(m.group(1).strip(), (m.group(2), m.group(3), m.group(4)))
    for lab, pre in KEEP.items():
        n, bp, pc = rows.get(lab, (None, None, None))
        d[pre+"_n"], d[pre+"_bp"], d[pre+"_pct"] = n, bp, pc
    return d

def p_rmod(p):
    t = rd(p)
    fams = re.findall(r"Round Time:.*?Elapsed Time\s*:\s*(\d+)\s+families discovered", t)
    cov = re.findall(r"([\d.]+)\s*%\s*of database", t)
    return {
        "repeatmodeler_version": g1(r"RepeatModeler Version\s+(\S+)", t),
        "rmod_search_engine":    g1(r"^Search Engine\s*=\s*(.+)$", t),
        "rmod_threads":          g1(r"^Threads\s*=\s*(\d+)", t),
        "rmod_dependencies":     g1(r"^Dependencies:\s*(.+)$", t),
        "ltr_module":            g1(r"LTR Structural Analysis:\s*(\w+)", t),
        "ltr_families":          g1(r"(\d+)\s+LTR families discovered", t),
        "n_rounds":              len(fams) or None,
        "families_per_round":    "/".join(fams) or None,
        "total_families_final":  g1(r"discovery complete:\s*(\d+)\s+families found", t),
        "final_coverage_pct":    cov[-1] if cov else None,
        "rmod_time":             g1(r"Program Time:\s*(\d+:\d+:\d+)", t),
        "rm_run_dir":            g1(r"(RM_\d+\S*)", t),
    }

def p_fam(p):
    cls = []
    for ln in open(p, "r", errors="replace"):
        if ln.startswith(">"):
            h = ln[1:].strip()
            cls.append(h.split("#", 1)[1].split()[0] if "#" in h else "NoClass")
    c = defaultdict(int)
    for x in cls: c[x] += 1
    tot = len(cls)
    unk = sum(v for k, v in c.items() if k.lower().startswith("unknown"))
    return {"consensi_total": tot, "consensi_unknown": unk,
            "consensi_classified": tot-unk,
            "consensi_unknown_frac": round(unk/tot, 4) if tot else None,
            "consensi_composition": "; ".join(f"{k}:{v}" for k, v in
                sorted(c.items(), key=lambda kv: -kv[1])[:12])}

def p_log(p):
    t = rd(p)
    return {"repeatmasker_version": g1(r"RepeatMasker\s+version\s+(\S+)", t)
                                    or g1(r"RepeatMasker\s+open-([\d.\-]+)", t),
            "softmask_xsmall": "Yes" if "-xsmall" in t else "No",
            "lib_arg": g1(r"-lib\s+(\S+)", t),
            "pipeline_wallclock": g1(r"(\d+:\d{2}:\d{2})(?!.*\d+:\d{2}:\d{2})", t, re.S)}

a = argparse.ArgumentParser()
for o in ("tbl-dir","rmod-dir","families-dir","log-dir"): a.add_argument("--"+o)
a.add_argument("--outdir", default="./S5_output")
A = a.parse_args()
P = lambda x: os.path.expanduser(x) if x else x
tb, rm, fa, lg = index(P(A.tbl_dir), [".tbl"]), index(P(A.rmod_dir), [".log"]), \
                 index(P(A.families_dir), [".fa", ".fasta"]), index(P(A.log_dir), [".log", ".txt"])
print(f".tbl {len(tb)} | rmod {len(rm)} | families {len(fa)} | log {len(lg)}")
os.makedirs(P(A.outdir), exist_ok=True)

recs, warn = [], []
for s in sorted(tb):
    r = OrderedDict(sample=s); miss = []
    r.update(p_tbl(tb[s]))
    for lab, idx, fn in (("rmod", rm, p_rmod), ("families", fa, p_fam), ("log", lg, p_log)):
        q = look(idx, s)
        if q:
            try: r.update(fn(q))
            except Exception as e: miss.append(f"{lab}:err({e})")
        else: miss.append(f"{lab}:missing")
    for k in ("total_masked_bp","interspersed_bp","consensi_total",
              "families_per_round","total_families_final","repeatmasker_version","rm_library"):
        if not r.get(k): miss.append("empty:"+k)
    r["parse_warnings"] = "; ".join(miss)
    recs.append(r)
    if miss: warn.append((s, "; ".join(miss)))

cols = []
for r in recs:
    for k in r:
        if k not in cols: cols.append(k)
o1 = os.path.join(P(A.outdir), "S5_master_table.tsv")
with open(o1, "w") as f:
    f.write("\t".join(cols)+"\n")
    for r in recs:
        f.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols)+"\n")
o2 = os.path.join(P(A.outdir), "S5_parse_report.txt")
with open(o2, "w") as f:
    f.write(f"total: {len(recs)}\nwarned: {len(warn)}\n\n")
    for s, m in warn: f.write(f"{s}\t{m}\n")
print(f"[written] {o1} ({len(recs)} x {len(cols)})")
print(f"[written] {o2} (warned {len(warn)})")
