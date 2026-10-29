#!/usr/bin/env python3
# 外源序列过滤后的功能注释重算（只读）
# 对 10 株过滤版蛋白集，只保留过滤后仍存在的蛋白 ID，重算功能注释统计，并替换进 202 株总表
# 用法见文末 README 或回复
import sys, os, glob, gzip, csv, re, collections, traceback, argparse
# [release refactor] Historical options: --freeze DIR (eggNOG/IPS = first glob match of
# emapper_v1_* / interproscan_v1_* inside it), --filtered, --table, --orig; output was
# dirname(--table)/func_filtered10. Directories are now explicit; statistics unchanged.
_ap = argparse.ArgumentParser(description="Recompute functional statistics for post-filter proteomes and merge into the 202-strain table.")
_ap.add_argument("--emapper-dir", required=True, help="frozen eggNOG directory (historical: <freeze>/emapper_v1_*)")
_ap.add_argument("--ips-dir", required=True, help="frozen InterProScan directory (historical: <freeze>/interproscan_v1_*)")
_ap.add_argument("--filtered", required=True, help="directory with post-filter <strain>_longest_isoform*.faa")
_ap.add_argument("--table", required=True, help="01_per_strain_functional_summary.tsv from summarize_functional_annotation.py (pre-filter)")
_ap.add_argument("--orig", default=None, help="optional: pre-filter <strain>_longest_isoform.faa directory for ID cross-check")
_ap.add_argument("--out", required=True, help="output directory (must not exist or be empty)")
_args = _ap.parse_args()
FD, TB, OR = _args.filtered, _args.table, _args.orig
for _d in [_args.emapper_dir, _args.ips_dir, FD] + ([OR] if OR else []):
    if not os.path.isdir(_d): sys.exit(f"[ERROR] directory not found: {_d}")
if not os.path.isfile(TB): sys.exit(f"[ERROR] table not found: {TB}")
OUT = _args.out
if os.path.isdir(OUT) and os.listdir(OUT): sys.exit(f"[ERROR] output directory not empty: {OUT}")
os.makedirs(OUT, exist_ok=True)
_STATUS = {"ok": False, "unmatched": 0}
log = open(os.path.join(OUT, "00_log.txt"), "w")
def say(*x):
    s = " ".join(map(str, x)); print(s); log.write(s+"\n"); log.flush()
key = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())
ids = lambda f: [l[1:].split()[0] for l in open(f) if l.startswith(">")]
APPL = ["Pfam","PANTHER","SUPERFAMILY","CDD","SMART","ProSitePatterns","ProSiteProfiles","AntiFam"]
def stats(s, keep, EM, IP):
    n=len(keep); ips=set(); per={x:set() for x in APPL}; ipr=set(); ig=set(); igi=set(); igp=set()
    with gzip.open(os.path.join(IP, s+".interproscan.tsv.gz"), "rt") as fh:
        for line in fh:
            c=line.rstrip("\n").split("\t")
            if len(c)<5 or c[0] not in keep: continue
            p=c[0]; ips.add(p)
            if c[3] in per: per[c[3]].add(p)
            if len(c)>11 and c[11] not in ("","-"): ipr.add(p)
            if len(c)>13 and c[13] not in ("","-"):
                ig.add(p)
                if "(InterPro)" in c[13]: igi.add(p)
                if "(PANTHER)" in c[13]: igp.add(p)
    egg=set(); e=collections.defaultdict(set); ix={}
    for line in open(os.path.join(EM, s+".emapper.annotations"), errors="replace"):
        if line.startswith("#query"): ix={k:j for j,k in enumerate(line[1:].rstrip("\n").split("\t"))}; continue
        if line.startswith("#") or not line.strip(): continue
        c=line.rstrip("\n").split("\t"); p=c[0]
        if p not in keep: continue
        egg.add(p); v=lambda k: c[ix[k]] if k in ix and ix[k]<len(c) else "-"
        for k,kk in [("GOs","go"),("KEGG_ko","ko"),("EC","ec"),("CAZy","cazy"),("KEGG_Pathway","path"),("Preferred_name","name"),("COG_category","cog")]:
            if v(k) not in ("-",""): e[kk].add(p)
    un=ips|egg; gu=ig|e["go"]; pc=lambda x: round(100*len(x)/n,2) if n else 0
    return dict(strain=s, n_proteins=n, ips_hit=len(ips), ips_pct=pc(ips), ips_ipr=len(ipr), ips_ipr_pct=pc(ipr),
        **{f"ips_{x}":len(per[x]) for x in APPL}, ips_go=len(ig), ips_go_pct=pc(ig), ips_go_InterPro=len(igi), ips_go_PANTHER=len(igp),
        egg_hit=len(egg), egg_pct=pc(egg), egg_go=len(e["go"]), egg_go_pct=pc(e["go"]), egg_ko=len(e["ko"]), egg_ko_pct=pc(e["ko"]),
        egg_ec=len(e["ec"]), egg_cazy=len(e["cazy"]), egg_cog=len(e["cog"]), egg_kegg_pathway=len(e["path"]), egg_preferred_name=len(e["name"]),
        union_hit=len(un), union_pct=pc(un), go_union=len(gu), go_union_pct=pc(gu),
        both_tools=len(ips&egg), ips_only=len(ips-egg), egg_only=len(egg-ips))
def main():
    EM=os.path.abspath(_args.emapper_dir); IP=os.path.abspath(_args.ips_dir)
    say("eggNOG:",EM); say("InterProScan:",IP); say("过滤版目录:",FD); say("原始 faa 目录:",OR)
    base=list(csv.DictReader(open(TB),delimiter="\t")); H=list(base[0].keys()); strains=[r["strain"] for r in base]
    fa=sorted(glob.glob(os.path.join(FD,"*.faa")))
    say("过滤版蛋白文件:",len(fa)); match={}
    for f in fa:
        b=os.path.basename(f); s=b.split("_longest_isoform")[0]
        if s in strains: match[s]=f
        else: say("  !! 未匹配:",b); _STATUS["unmatched"]+=1
    new={}; cmp=[]
    for s,f in sorted(match.items()):
        kl=ids(f); keep=set(kl)
        if len(kl)!=len(keep): say(f"  !! {s} 过滤版中有重复 ID")
        if OR:
            o=os.path.join(OR, s+"_longest_isoform.faa")
            if os.path.exists(o):
                orig=set(ids(o)); extra=keep-orig
                say(f"  {s}: 原始 {len(orig)}，过滤版 {len(keep)}，移除 {len(orig-keep)}，过滤版中不在原始集的 ID {len(extra)}")
                if extra: say("  !! 有 ID 不在原始集中，示例:", list(extra)[:3])
        b0=next(r for r in base if r["strain"]==s); a1=stats(s,keep,EM,IP); new[s]=a1
        cmp.append([s,b0["n_proteins"],a1["n_proteins"],int(b0["n_proteins"])-a1["n_proteins"],b0["ips_pct"],a1["ips_pct"],b0["egg_pct"],a1["egg_pct"],
                    b0["union_pct"],a1["union_pct"],b0["go_union_pct"],a1["go_union_pct"],b0["ips_AntiFam"],a1["ips_AntiFam"]])
        say(f"    IPS {b0['ips_pct']}→{a1['ips_pct']}  eggNOG {b0['egg_pct']}→{a1['egg_pct']}  合计 {b0['union_pct']}→{a1['union_pct']}  AntiFam {b0['ips_AntiFam']}→{a1['ips_AntiFam']}")
    with open(os.path.join(OUT,"01_filtered10_before_after.tsv"),"w") as o:
        o.write("strain\tn_before\tn_after\tn_removed\tips_pct_before\tips_pct_after\tegg_pct_before\tegg_pct_after\tunion_pct_before\tunion_pct_after\tgo_union_pct_before\tgo_union_pct_after\tantifam_before\tantifam_after\n")
        for r in cmp: o.write("\t".join(map(str,r))+"\n")
    fin=[new.get(r["strain"],r) for r in base]
    with open(os.path.join(OUT,"02_per_strain_functional_summary_final.tsv"),"w") as o:
        o.write("\t".join(H+["filtered"])+"\n")
        for r in fin: o.write("\t".join(str(r[k]) for k in H)+"\t"+("yes" if r["strain"] in new else "no")+"\n")
    T=sum(int(r["n_proteins"]) for r in fin); say("== 最终版本（202 株，10 株为过滤版）蛋白总数:",T)
    for k in ["ips_hit","ips_ipr","ips_go","ips_go_InterPro","ips_go_PANTHER","egg_hit","egg_go","egg_ko","egg_ec","egg_cazy","egg_cog","egg_kegg_pathway","union_hit","go_union","ips_AntiFam","both_tools","ips_only","egg_only"]:
        t=sum(int(r[k]) for r in fin); say(f"  {k}: {t}（{100*t/T:.2f}%）")
    say(f"  AntiFam 有命中的株: {sum(1 for r in fin if int(r['ips_AntiFam'])>0)}")
    say("完成。请把 func_filtered10 文件夹传给我")
    _STATUS["ok"] = True
try: main()
except Exception: say("!! 出错:\n"+traceback.format_exc())
# [release refactor] historical script always exited 0
if not _STATUS["ok"]: sys.exit(1)
if _STATUS["unmatched"]: say(f"!! {_STATUS['unmatched']} 个过滤版文件未匹配到总表株名"); sys.exit(2)
