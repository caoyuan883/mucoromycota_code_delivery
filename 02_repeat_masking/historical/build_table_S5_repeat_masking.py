#!/usr/bin/env python3
"""
build_table_S5_repeat_masking.py

汇总 RepeatModeler / RepeatMasker 输出, 生成表5 主表 (每株一行) 与可选的
重复景观长表 (genome x class x divergence-bin -> bp)。

依赖: 仅 Python 3.8+ 标准库 (pandas 可选, 无则输出 TSV 仍正常)。

用法:
  python3 build_table_S5_repeat_masking.py \
      --tbl-dir      ~/Desktop/211_mucoromycota_repeat_and_logs/05_collected_key_outputs/repeatmasker_tbl \
      --rmod-dir     ~/Desktop/211_mucoromycota_repeat_and_logs/05_collected_key_outputs/rmod_logs \
      --families-dir ~/Desktop/211_mucoromycota_repeat_and_logs/05_collected_key_outputs/families_fa \
      --log-dir      ~/Desktop/211_mucoromycota_repeat_and_logs/05_parallel_logs \
      --out-dir      /path/to/dir_with_fna_out            # 可选, 有则生成景观表
      --outdir       ./S5_output

先跑一次体检 (强烈建议):
  python3 build_table_S5_repeat_masking.py --inspect --rmod-dir <...> --log-dir <...>
它会打印 rmod.log / parallel log 里的关键行, 便于你确认下方 PATTERNS 是否匹配你的版本。
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict, OrderedDict

# ---------------------------------------------------------------- 可调正则 --
# 不同 RepeatModeler 版本措辞略有差异; 如 --inspect 显示不匹配, 在此增补即可。
PATTERNS = {
    "rmodeler_version": [r"RepeatModeler\s+Version\s+(\S+)", r"RepeatModeler\s+(\d+\.\d+\.\d+)"],
    "search_engine":    [r"Search Engine\s*=\s*(.+)", r"Using search engine\s*[:=]\s*(.+)"],
    "threads":          [r"-pa\s+(\d+)", r"Number of threads\s*[:=]\s*(\d+)", r"threads\s*[:=]\s*(\d+)"],
    "round_families":   [r"-\s*(\d+)\s+families\s+discovered",
                         r"Number of families returned by RepeatScout:\s*(\d+)",
                         r"-\s*Ended with\s+(\d+)\s+families"],
    "round_marker":     [r"^\s*Round\s*#?\s*(\d+)"],
    "coverage":         [r"\(\s*([\d.]+)\s*%\s*of database\s*\)", r"[Cc]overage\s*[:=]\s*([\d.]+)\s*%"],
    "rmod_time":        [r"Program Time:\s*([\d:]+)", r"Total time\s*[:=]\s*([\d:]+)"],
    "rm_dir":           [r"(RM_\d+\.\S*)", r"(RM_\d+)"],
    "rm_version":       [r"RepeatMasker\s+version\s+(\S+)", r"RepeatMasker\s+open-([\d.\-]+)"],
    "xsmall":           [r"(-xsmall)"],
    "lib_arg":          [r"-lib\s+(\S+)"],
    "wallclock":        [r"real\s+(\S+)", r"[Ee]lapsed[^\d]*([\d:]+)", r"Total run time\s*[:=]\s*([\d:]+)"],
}

TBL_ROW_RE = re.compile(r"^\s*([A-Za-z][^:]*?):?\s+(\d+)\s+(\d+)\s*bp\s+([\d.]+)\s*%\s*$")

# .tbl 中希望固定输出的行标签 -> 输出列前缀
TBL_KEEP = OrderedDict([
    ("Retroelements",      "retroelement"),
    ("SINEs",              "SINE"),
    ("LINEs",              "LINE"),
    ("LTR elements",       "LTR"),
    ("DNA transposons",    "DNA_TE"),
    ("Rolling-circles",    "RC_Helitron"),
    ("Unclassified",       "unclassified"),
    ("Small RNA",          "small_RNA"),
    ("Satellites",         "satellite"),
    ("Simple repeats",     "simple_repeat"),
    ("Low complexity",     "low_complexity"),
])


def sopen(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") \
        else open(path, "r", errors="replace")


def first_match(text, keys):
    for pat in PATTERNS[keys]:
        m = re.search(pat, text, flags=re.M)
        if m:
            return m.group(1).strip()
    return None


def norm_sample(fname):
    """从文件名剥离常见后缀, 得到样本名。"""
    b = os.path.basename(fname)
    for suf in (".gz",):
        if b.endswith(suf):
            b = b[: -len(suf)]
    for suf in (".fna.tbl", ".fa.tbl", ".tbl", "-rmod.log", "_rmod.log",
                "-families.fa", "_families.fa", ".fna.out", ".out", ".log"):
        if b.endswith(suf):
            b = b[: -len(suf)]
            break
    b = re.sub(r"_genomic$", "", b)
    b = re.sub(r"[._-]+$", "", b)
    return b


def index_dir(d, exts):
    """返回 {sample_name: path}"""
    out = {}
    if not d or not os.path.isdir(d):
        return out
    for fn in sorted(os.listdir(d)):
        if any(fn.endswith(e) or fn.endswith(e + ".gz") for e in exts):
            out.setdefault(norm_sample(fn), os.path.join(d, fn))
    return out


def lookup(index, sample):
    """精确 -> 前缀 -> 去分隔符 模糊匹配。"""
    if sample in index:
        return index[sample]
    for k, v in index.items():
        if k.startswith(sample) or sample.startswith(k):
            return v
    flat = re.sub(r"[^a-z0-9]", "", sample.lower())
    for k, v in index.items():
        if re.sub(r"[^a-z0-9]", "", k.lower()) == flat:
            return v
    return None


# ------------------------------------------------------------------ 解析器 --
def parse_tbl(path):
    with sopen(path) as fh:
        txt = fh.read()
    d = {}
    m = re.search(r"sequences:\s+(\d+)", txt)
    d["sequences"] = int(m.group(1)) if m else None
    m = re.search(r"total length:\s+(\d+)\s*bp", txt)
    d["total_length_bp"] = int(m.group(1)) if m else None
    m = re.search(r"GC level:\s+([\d.]+)\s*%", txt)
    d["GC_percent"] = float(m.group(1)) if m else None
    m = re.search(r"bases masked:\s+(\d+)\s*bp\s*\(\s*([\d.]+)\s*%\)", txt)
    d["total_masked_bp"] = int(m.group(1)) if m else None
    d["total_masked_pct"] = float(m.group(2)) if m else None
    m = re.search(r"Total interspersed repeats:\s+(\d+)\s*bp\s+([\d.]+)\s*%", txt)
    d["interspersed_bp"] = int(m.group(1)) if m else None
    d["interspersed_pct"] = float(m.group(2)) if m else None
    m = re.search(r"RepeatMasker Combined Database:\s*(\S+)", txt)
    d["rm_library"] = m.group(1) if m else None
    m = re.search(r"run with\s+(\S+\s+version\s+\S+)", txt)
    d["rm_search_engine"] = m.group(1) if m else None

    rows = {}
    for line in txt.splitlines():
        mm = TBL_ROW_RE.match(line)
        if mm:
            rows.setdefault(mm.group(1).strip(),
                            (int(mm.group(2)), int(mm.group(3)), float(mm.group(4))))
    for label, pref in TBL_KEEP.items():
        n, bp, pct = rows.get(label, (None, None, None))
        d[pref + "_n"] = n
        d[pref + "_bp"] = bp
        d[pref + "_pct"] = pct
    d["_tbl_all_rows"] = rows
    return d


def parse_rmod(path):
    with sopen(path) as fh:
        txt = fh.read()
    d = {
        "repeatmodeler_version": first_match(txt, "rmodeler_version"),
        "rmod_search_engine":    first_match(txt, "search_engine"),
        "rmod_threads":          first_match(txt, "threads"),
        "rmod_time":             first_match(txt, "rmod_time"),
    }
    fams = []
    for pat in PATTERNS["round_families"]:
        fams = re.findall(pat, txt, flags=re.M)
        if fams:
            break
    d["families_per_round"] = "/".join(fams) if fams else None
    d["n_rounds"] = len(fams) if fams else None
    d["total_families_log"] = sum(int(x) for x in fams) if fams else None
    cov = []
    for pat in PATTERNS["coverage"]:
        cov = re.findall(pat, txt, flags=re.M)
        if cov:
            break
    d["final_coverage_pct"] = cov[-1] if cov else None
    m = re.search(PATTERNS["rm_dir"][0], txt)
    d["rm_run_dir"] = m.group(1) if m else None
    d["dfam_warning"] = "Yes" if re.search(r"Dfam.*(incomplete|not.*complete|partial)",
                                           txt, flags=re.I) else "No"
    return d


def parse_families(path):
    with sopen(path) as fh:
        heads = [ln[1:].strip() for ln in fh if ln.startswith(">")]
    classes = []
    for h in heads:
        classes.append(h.split("#", 1)[1].split()[0] if "#" in h else "NoClass")
    cnt = defaultdict(int)
    for c in classes:
        cnt[c] += 1
    total = len(classes)
    unknown = sum(v for k, v in cnt.items() if k.lower().startswith("unknown"))
    top = "; ".join(f"{k}:{v}" for k, v in
                    sorted(cnt.items(), key=lambda kv: -kv[1])[:12])
    return {
        "consensi_total": total,
        "consensi_unknown": unknown,
        "consensi_unknown_frac": round(unknown / total, 4) if total else None,
        "consensi_classified": total - unknown,
        "consensi_composition": top,
    }


def parse_pipeline_log(path):
    with sopen(path) as fh:
        txt = fh.read()
    return {
        "repeatmasker_version": first_match(txt, "rm_version"),
        "softmask_xsmall":      "Yes" if first_match(txt, "xsmall") else "No/未检出",
        "lib_arg":              first_match(txt, "lib_arg"),
        "rm_run_dir_log":       first_match(txt, "rm_dir"),
        "pipeline_wallclock":   first_match(txt, "wallclock"),
    }


def parse_out(path, bin_size=1.0, max_div=60.0):
    """返回 (per-class 汇总, 景观 bins)。divergence 取 .out 第 2 列 perc div.。"""
    per_class_bp = defaultdict(int)
    per_class_wdiv = defaultdict(float)
    bins = defaultdict(int)          # (class, bin_start) -> bp
    n = 0
    with sopen(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) < 14 or not t[0].isdigit():
                continue
            try:
                div = float(t[1]); qb = int(t[5]); qe = int(t[6]); cls = t[10]
            except (ValueError, IndexError):
                continue
            ln = qe - qb + 1
            if ln <= 0:
                continue
            n += 1
            per_class_bp[cls] += ln
            per_class_wdiv[cls] += div * ln
            b = min(div, max_div)
            bins[(cls, round((b // bin_size) * bin_size, 1))] += ln
    summ = {"out_n_hits": n}
    tot_bp = sum(per_class_bp.values())
    if tot_bp:
        summ["out_total_hit_bp"] = tot_bp
        summ["out_mean_div_bp_weighted"] = round(
            sum(per_class_wdiv.values()) / tot_bp, 3)
        for major in ("LTR", "LINE", "DNA", "Unknown", "RC"):
            bp = sum(v for k, v in per_class_bp.items() if k.startswith(major))
            wd = sum(per_class_wdiv[k] for k in per_class_bp if k.startswith(major))
            summ[f"div_mean_{major}"] = round(wd / bp, 3) if bp else None
    return summ, bins


# --------------------------------------------------------------------- 主程 --
def inspect(args):
    for label, d, exts in (("rmod.log", args.rmod_dir, [".log"]),
                           ("pipeline log", args.log_dir, [".log", ".txt"])):
        idx = index_dir(d, exts)
        if not idx:
            print(f"[{label}] 目录为空或不存在: {d}")
            continue
        k = sorted(idx)[0]
        print(f"\n===== [{label}] 样例文件: {idx[k]} =====")
        with sopen(idx[k]) as fh:
            for ln in fh:
                if re.search(r"Version|Round|famil|Coverage|of database|Time|"
                             r"RM_\d+|xsmall|-lib|real\s|Engine|thread",
                             ln, flags=re.I):
                    print("   " + ln.rstrip()[:160])
    print("\n若上面的关键行与脚本顶部 PATTERNS 不匹配, 直接在 PATTERNS 里增补正则即可。")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tbl-dir", required=False)
    p.add_argument("--rmod-dir", required=False)
    p.add_argument("--families-dir", required=False)
    p.add_argument("--log-dir", required=False)
    p.add_argument("--out-dir", required=False, help="含 .fna.out 的目录 (可选)")
    p.add_argument("--outdir", default="./S5_output")
    p.add_argument("--bin-size", type=float, default=1.0)
    p.add_argument("--inspect", action="store_true")
    a = p.parse_args()

    for k in ("tbl_dir", "rmod_dir", "families_dir", "log_dir", "out_dir", "outdir"):
        v = getattr(a, k)
        if v:
            setattr(a, k, os.path.expanduser(v))

    if a.inspect:
        inspect(a); return

    if not a.tbl_dir:
        sys.exit("必须提供 --tbl-dir (样本清单以 .tbl 为准)")

    tbl_idx = index_dir(a.tbl_dir, [".tbl"])
    rmod_idx = index_dir(a.rmod_dir, [".log"])
    fam_idx = index_dir(a.families_dir, [".fa", ".fasta"])
    log_idx = index_dir(a.log_dir, [".log", ".txt"])
    out_idx = index_dir(a.out_dir, [".out"])

    print(f"发现 .tbl {len(tbl_idx)} 份 | rmod.log {len(rmod_idx)} | "
          f"families.fa {len(fam_idx)} | pipeline log {len(log_idx)} | .out {len(out_idx)}")

    os.makedirs(a.outdir, exist_ok=True)
    records, missing, landscape_rows = [], [], []

    for sample in sorted(tbl_idx):
        rec = OrderedDict(sample=sample)
        miss = []

        d = parse_tbl(tbl_idx[sample]); d.pop("_tbl_all_rows", None)
        rec.update(d)

        for label, idx, fn in (("rmod.log", rmod_idx, parse_rmod),
                               ("families.fa", fam_idx, parse_families),
                               ("pipeline.log", log_idx, parse_pipeline_log)):
            path = lookup(idx, sample)
            if path:
                try:
                    rec.update(fn(path))
                except Exception as ex:            # noqa: BLE001
                    miss.append(f"{label}:解析失败({ex})")
            else:
                miss.append(f"{label}:缺失")

        if out_idx:
            path = lookup(out_idx, sample)
            if path:
                summ, bins = parse_out(path, a.bin_size)
                rec.update(summ)
                for (cls, b), bp in bins.items():
                    landscape_rows.append((sample, cls, b, bp))
            else:
                miss.append(".out:缺失")

        for key in ("total_masked_bp", "interspersed_bp", "consensi_total",
                    "families_per_round", "repeatmasker_version"):
            if rec.get(key) in (None, ""):
                miss.append(f"字段缺失:{key}")

        rec["parse_warnings"] = "; ".join(miss) if miss else ""
        records.append(rec)
        if miss:
            missing.append((sample, "; ".join(miss)))

    cols = []
    for r in records:
        for k in r:
            if k not in cols:
                cols.append(k)

    master = os.path.join(a.outdir, "S5_master_table.tsv")
    with open(master, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c))
                               for c in cols) + "\n")
    print(f"[写出] {master}  ({len(records)} 行 x {len(cols)} 列)")

    if landscape_rows:
        ls = os.path.join(a.outdir, "S5_repeat_landscape_long.tsv")
        with open(ls, "w") as fh:
            fh.write("sample\trepeat_class\tdiv_bin_start\tbp\n")
            for row in landscape_rows:
                fh.write("\t".join(map(str, row)) + "\n")
        print(f"[写出] {ls}  ({len(landscape_rows)} 行)")

    rep = os.path.join(a.outdir, "S5_parse_report.txt")
    with open(rep, "w") as fh:
        fh.write(f"样本总数: {len(records)}\n有告警的样本: {len(missing)}\n\n")
        for s, m in missing:
            fh.write(f"{s}\t{m}\n")
    print(f"[写出] {rep}  (有告警样本 {len(missing)} 个 —— 请务必先看这个文件)")


if __name__ == "__main__":
    main()
