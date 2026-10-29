#!/usr/bin/env python3
# 功能注释汇总 v3：直接从 func_freeze 的结果文件和日志统计；自动跳过缺失、MD5 不一致和损坏的文件
# 用法 (release): python3 summarize_functional_annotation.py --emapper-dir DIR --ips-dir DIR --faa DIR --out DIR
# [release refactor] Historical usage was `func_summary_v3.py [base] --faa DIR` with base
# defaulting to ~/Desktop/interproscan_eggnog; the eggNOG/IPS directories were the sorted-LAST
# match of emapper_v1_* / interproscan_v1_* under base, and without --faa the proteins were
# searched under base and then ~/Desktop (first match won). These implicit choices are now
# explicit arguments. Statistics (including the nearest-rank median) are unchanged.
import sys, os, glob, gzip, hashlib, csv, collections, time, re, traceback, argparse
_ap = argparse.ArgumentParser(description="Per-strain InterProScan/eggNOG functional summary (frozen result directories).")
_ap.add_argument("--emapper-dir", required=True, help="frozen eggNOG directory (historical: emapper_v1_*)")
_ap.add_argument("--ips-dir", required=True, help="frozen InterProScan directory (historical: interproscan_v1_*)")
_ap.add_argument("--faa", required=True, help="directory with <strain>_longest_isoform.faa (defines n_proteins denominators)")
_ap.add_argument("--out", required=True, help="output directory (must not exist or be empty)")
_args = _ap.parse_args()
FAA = _args.faa
for _d in (_args.emapper_dir, _args.ips_dir, FAA):
    if not os.path.isdir(_d): sys.exit(f"[ERROR] directory not found: {_d}")
OUT = _args.out
if os.path.isdir(OUT) and os.listdir(OUT): sys.exit(f"[ERROR] output directory not empty: {OUT}")
os.makedirs(OUT, exist_ok=True)
_STATUS = {"ok": False}
log = open(os.path.join(OUT, "00_run_log.txt"), "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()
def wtsv(name, header, data):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(data)
def files(d, pat, pre, suf):
    out = {}
    for f in glob.glob(os.path.join(d, "**", pat), recursive=True):
        b = os.path.basename(f)
        if ".stuck_" in b: continue
        out[b[len(pre):len(b)-len(suf)]] = f
    return out
def sec(t):
    if not t: return ""
    p = [float(x) for x in t.split(":")]
    return int(p[0]*3600 + p[1]*60 + p[2]) if len(p) == 3 else int(p[0]*60 + p[1])
def last(txt, pat):
    m = re.findall(pat, txt); return m[-1] if m else ""
ELAP = r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)"
APPL = ["Pfam","PANTHER","SUPERFAMILY","CDD","SMART","ProSitePatterns","ProSiteProfiles","AntiFam"]

def main():
    EM, IP = os.path.abspath(_args.emapper_dir), os.path.abspath(_args.ips_dir)
    say("eggNOG 目录:", EM); say("InterProScan 目录:", IP); say("输出目录:", OUT)
    if not EM or not IP: say("!! 没找到两个定版目录"); return
    em_ann = files(EM, "*.emapper.annotations", "", ".emapper.annotations")
    ip_tsv = files(IP, "*.interproscan.tsv.gz", "", ".interproscan.tsv.gz")
    em_log = files(EM, "em_*.log", "em_", ".log"); ip_log = files(IP, "ips_*.log", "ips_", ".log")
    say(f"== 0. 现有文件: eggNOG 注释 {len(em_ann)}，IPS TSV {len(ip_tsv)}，eggNOG 日志 {len(em_log)}，IPS 日志 {len(ip_log)}")

    say("== 1. MD5 校验")
    bad = set()
    for d in (EM, IP):
        f = os.path.join(d, "MD5SUMS")
        if not os.path.exists(f): say("  无 MD5SUMS:", os.path.basename(d)); continue
        n = miss = mis = 0
        for line in open(f):
            if not line.strip(): continue
            h, fn = line.strip().split(None, 1); fn = fn.lstrip("*"); p = os.path.join(d, fn); n += 1
            if not os.path.exists(p): miss += 1; continue
            m = hashlib.md5()
            with open(p, "rb") as fh:
                for ch in iter(lambda: fh.read(1 << 22), b""): m.update(ch)
            if m.hexdigest() != h: mis += 1; bad.add(os.path.abspath(p)); say("  MD5 不一致（跳过）:", fn)
        say(f"  {os.path.basename(d)}: 清单 {n}，缺失 {miss}，不一致 {mis}")
    for dct in (em_ann, ip_tsv):
        for k in [k for k, f in dct.items() if os.path.abspath(f) in bad]: del dct[k]

    say("== 2. 蛋白数（来自 faa）")
    nprot = {}
    for c in [FAA]:
        if not c or not os.path.isdir(c): continue
        for f in glob.glob(os.path.join(c, "**", "*_longest_isoform.faa"), recursive=True):
            b = os.path.basename(f)[:-len("_longest_isoform.faa")]
            if b not in nprot: nprot[b] = sum(1 for l in open(f) if l.startswith(">"))
        if nprot: say(f"  来源 {c}，{len(nprot)} 株"); break
    strains = sorted(set(nprot) | set(em_log) | set(ip_log) | set(em_ann) | set(ip_tsv))
    say(f"  株总数 {len(strains)}；缺蛋白数 {sum(1 for s in strains if s not in nprot)}")

    say("== 3. 日志（不受下载缺失影响）")
    H6 = ["strain","n_proteins","ips_unique_seqs","ips_100pct_done","ips_exit","ips_wall_s","ips_maxrss_kb","ips_cpu_pct",
          "ips_warn_failed_match","ips_warn_npe","egg_hits_processed","egg_pct_from_log","egg_finished","egg_exit",
          "egg_wall_s","egg_maxrss_kb","egg_cpu_pct","ips_log"]
    lr = []
    for s in strains:
        ii = open(ip_log[s], errors="replace").read() if s in ip_log else ""
        ei = open(em_log[s], errors="replace").read() if s in em_log else ""
        n = nprot.get(s); hp = last(ei, r"Total hits processed: (\d+)")
        lr.append([s, n or "NA", last(ii, r"Uploaded (\d+) unique sequences"),
                   "yes" if "100% done" in ii else ("no" if ii else "NA"), last(ii, r"Exit status: (\d+)"),
                   sec(last(ii, ELAP)), last(ii, r"Maximum resident set size \(kbytes\): (\d+)"), last(ii, r"Percent of CPU this job got: (\d+)%"),
                   ii.count("Failed to get match results"), ii.count("NullPointerException"),
                   hp, round(100*int(hp)/n, 2) if hp and n else "NA",
                   "yes" if "FINISHED" in ei else ("no" if ei else "NA"), last(ei, r"Exit status: (\d+)"),
                   sec(last(ei, ELAP)), last(ei, r"Maximum resident set size \(kbytes\): (\d+)"), last(ei, r"Percent of CPU this job got: (\d+)%"),
                   "normal" if ii else "none(rerun?)"])
    wtsv("06_run_logs_per_strain.tsv", H6, lr)
    say(f"  IPS 日志 {sum(r[3]!='NA' for r in lr)} 株：100% done {sum(r[3]=='yes' for r in lr)}，退出码0 {sum(r[4]=='0' for r in lr)}，丢失警告>0 {sum(r[8]>0 for r in lr)} 株，NPE>0 {sum(r[9]>0 for r in lr)} 株")
    for r in lr:
        if r[8] > 0 or r[9] > 0: say(f"    警告 {r[0]}: failed_match {r[8]}，NPE {r[9]}")
    say(f"  eggNOG 日志 {sum(r[12]!='NA' for r in lr)} 株：FINISHED {sum(r[12]=='yes' for r in lr)}，退出码0 {sum(r[13]=='0' for r in lr)}")
    ep = sorted(r[11] for r in lr if r[11] != "NA")
    if ep:
        q = lambda p: ep[min(len(ep)-1, int(round(p*(len(ep)-1))))]
        tot_h = sum(int(r[10]) for r in lr if r[10]); tot_n = sum(r[1] for r in lr if r[10] and r[1] != "NA")
        say(f"  eggNOG 覆盖率（日志 hits/蛋白数，{len(ep)} 株）: 合计 {100*tot_h/tot_n:.2f}%，中位数 {q(.5)}%（IQR {q(.25)}–{q(.75)}%），范围 {ep[0]}–{ep[-1]}%")
    wt = [r[5] for r in lr if r[5] != ""]; we = [r[14] for r in lr if r[14] != ""]
    if wt: say(f"  IPS 单株耗时: 中位数 {sorted(wt)[len(wt)//2]/3600:.1f} h，合计 {sum(wt)/3600:.0f} h（{len(wt)} 株）")
    if we: say(f"  eggNOG 单株耗时: 中位数 {sorted(we)[len(we)//2]/60:.0f} min，合计 {sum(we)/3600:.0f} h（{len(we)} 株）")

    say("== 4. 结果文件统计（只统计两个文件都完整的株）")
    rows = []; lvl = collections.Counter(); lvl_s = collections.defaultdict(collections.Counter)
    cog_s = collections.defaultdict(collections.Counter); gosrc = collections.Counter(); avail = []
    for s in strains:
        ok_e, ok_i = s in em_ann, s in ip_tsv
        avail.append([s, "yes" if ok_e else "no", "yes" if ok_i else "no"])
        if not (ok_e and ok_i and s in nprot): continue
        n = nprot[s]
        ips = set(); per = {a: set() for a in APPL}; ipr = set(); ig = set(); igi = set(); igp = set(); gs = collections.Counter()
        try:
            with gzip.open(ip_tsv[s], "rt") as fh:
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 5: continue
                    p = c[0]; ips.add(p)
                    if c[3] in per: per[c[3]].add(p)
                    if len(c) > 11 and c[11] not in ("", "-"): ipr.add(p)
                    if len(c) > 13 and c[13] not in ("", "-"):
                        ig.add(p)
                        for g in c[13].split("|"):
                            sr = g[g.find("(")+1:g.find(")")] if "(" in g else "unspecified"; gs[sr] += 1
                            if sr == "InterPro": igi.add(p)
                            elif sr == "PANTHER": igp.add(p)
        except (EOFError, OSError, gzip.BadGzipFile) as ex:
            say(f"  !! {s} IPS 文件损坏，跳过: {ex}"); avail[-1][2] = "corrupt"; continue
        egg = set(); e = collections.defaultdict(set); ix = {}; ls = collections.Counter(); cs = collections.Counter()
        for line in open(em_ann[s], errors="replace"):
            if line.startswith("#query"): ix = {k: j for j, k in enumerate(line[1:].rstrip("\n").split("\t"))}; continue
            if line.startswith("#") or not line.strip(): continue
            c = line.rstrip("\n").split("\t"); p = c[0]; egg.add(p)
            v = lambda k: c[ix[k]] if k in ix and ix[k] < len(c) else "-"
            for k, key in [("GOs","go"),("KEGG_ko","ko"),("EC","ec"),("CAZy","cazy"),("KEGG_Pathway","path"),("Preferred_name","name")]:
                if v(k) not in ("-", ""): e[key].add(p)
            cg = v("COG_category")
            if cg not in ("-", ""):
                e["cog"].add(p)
                for ch in cg: cs[ch] += 1
            ml = v("max_annot_lvl")
            if ml not in ("-", ""): ls[ml.split("|")[-1]] += 1
        gosrc.update(gs); lvl.update(ls); lvl_s[s] = ls; cog_s[s] = cs
        un = ips | egg; gu = ig | e["go"]; pc = lambda x: round(100*len(x)/n, 2)
        rows.append(dict(strain=s, n_proteins=n, ips_hit=len(ips), ips_pct=pc(ips), ips_ipr=len(ipr), ips_ipr_pct=pc(ipr),
            **{f"ips_{a}": len(per[a]) for a in APPL}, ips_go=len(ig), ips_go_pct=pc(ig), ips_go_InterPro=len(igi), ips_go_PANTHER=len(igp),
            egg_hit=len(egg), egg_pct=pc(egg), egg_go=len(e["go"]), egg_go_pct=pc(e["go"]), egg_ko=len(e["ko"]), egg_ko_pct=pc(e["ko"]),
            egg_ec=len(e["ec"]), egg_cazy=len(e["cazy"]), egg_cog=len(e["cog"]), egg_kegg_pathway=len(e["path"]), egg_preferred_name=len(e["name"]),
            union_hit=len(un), union_pct=pc(un), go_union=len(gu), go_union_pct=pc(gu),
            both_tools=len(ips & egg), ips_only=len(ips - egg), egg_only=len(egg - ips)))
        say(f"  {len(rows)} {s}  IPS {pc(ips)}%  eggNOG {pc(egg)}%  合计 {pc(un)}%")
    wtsv("07_file_availability.tsv", ["strain","egg_annotations","ips_tsv"], avail)
    say(f"  可统计 {len(rows)} 株（共 {len(strains)} 株）")
    if not rows: return
    H = list(rows[0]); wtsv("01_per_strain_functional_summary.tsv", H, [[r[k] for k in H] for r in rows])
    S = [r["strain"] for r in rows]; L = [k for k, _ in lvl.most_common()]
    wtsv("02_eggnog_max_annot_lvl_per_strain.tsv", ["strain"] + L, [[s] + [lvl_s[s].get(k, 0) for k in L] for s in S])
    C = sorted({ch for s in S for ch in cog_s[s]})
    wtsv("03_eggnog_COG_category_per_strain.tsv", ["strain"] + C, [[s] + [cog_s[s].get(ch, 0) for ch in C] for s in S])
    wtsv("04_ips_GO_source_counts.tsv", ["source","n_term_assignments"], gosrc.most_common())
    T = sum(r["n_proteins"] for r in rows)
    def st(k):
        v = sorted(r[k] for r in rows); q = lambda p: v[min(len(v)-1, int(round(p*(len(v)-1))))]
        return [k, q(.5), q(.25), q(.75), v[0], v[-1]]
    summ = [["strains", len(rows), "", "", "", ""], ["total_proteins", T, "", "", "", ""]]
    for k in ["ips_hit","ips_ipr","ips_go","egg_hit","egg_go","egg_ko","egg_ec","egg_cazy","egg_cog","union_hit","go_union","ips_AntiFam"]:
        t = sum(r[k] for r in rows); summ.append([k+"_total", t, round(100*t/T, 2), "", "", ""])
    for k in ["ips_pct","ips_ipr_pct","ips_go_pct","egg_pct","egg_go_pct","egg_ko_pct","union_pct","go_union_pct"]: summ.append(st(k))
    wtsv("05_cohort_summary.tsv", ["metric","median_or_total","q1_or_pct","q3","min","max"], summ)
    say("== 5. 汇总（仅可统计株）"); [say("  " + "\t".join(map(str, x))) for x in summ]
    say("  AntiFam 有命中的株:", sum(1 for r in rows if r["ips_AntiFam"] > 0))
    say("  eggNOG 注释层级前 8:", lvl.most_common(8))
    say("== 完成。请把 func_summary_v3 文件夹传给我")
    _STATUS["ok"] = True

try: main()
except Exception:
    say("!! 运行出错："); say(traceback.format_exc())
# [release refactor] historical script always exited 0; now non-zero unless all steps completed
if not _STATUS["ok"]:
    say("!! 未完成（缺少目录、无可统计株或运行出错）"); sys.exit(1)
