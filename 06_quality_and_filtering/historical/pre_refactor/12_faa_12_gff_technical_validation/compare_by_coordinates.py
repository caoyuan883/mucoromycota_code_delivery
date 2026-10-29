#!/usr/bin/env python3
import re, gzip, hashlib, statistics as st
from pathlib import Path
from collections import defaultdict

BASE = Path("/home/py/Desktop/12_faa_12_gff_technical_validation")
REF_FA, OWN_FA = BASE/"12_faa_longest_isoform", BASE/"12_faa_own_longest_isoform"
GFF, GTF = BASE/"12_gff", BASE/"12_gtf_own"
OUT = BASE/"coord_comparison"; OUT.mkdir(exist_ok=True)

MIN_OVL_BP, MIN_OVL_FRAC = 30, 0.20

def op(p): return gzip.open(p,"rt") if str(p).endswith(".gz") else open(p,"r")
def attrs(s):
    s=s.strip()
    if "=" in s and not re.search(r'\w+\s+"',s):
        d={}
        for kv in s.rstrip(";").split(";"):
            kv=kv.strip()
            if "=" in kv: k,v=kv.split("=",1); d[k.strip()]=v.strip()
        return d,"gff3"
    return {m.group(1):m.group(2) for m in re.finditer(r'(\S+)\s+"([^"]*)"',s)},"gtf"
def cands(v):
    if not v: return []
    v=v.strip(); out=[v,re.sub(r"\.\d+$","",v)]
    for pre in ("cds-","rna-","gene-"):
        if v.startswith(pre): t=v[len(pre):]; out+=[t,re.sub(r"\.\d+$","",t)]
    if "|" in v:
        f=[x for x in v.split("|") if x]
        for c in ([f[-1]] if f else [])+[x for x in f if x.isdigit()]:
            out+=[c,re.sub(r"\.\d+$","",c)]
    seen=set(); r=[]
    for x in out:
        if x and x not in seen: seen.add(x); r.append(x)
    return r

def load_ann(path):
    """-> gene_id: dict(seqid,strand,iv,cdslen,nex,keys)"""
    lines=[]
    with op(path) as fh:
        for line in fh:
            if not line or line[0]=="#": continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9: continue
            lines.append(f)
    m2g={}
    for f in lines:
        a,fmt=attrs(f[8])
        if fmt=="gff3" and f[2] in ("mRNA","transcript"):
            i=a.get("ID","")
            if i: m2g[i]=a.get("Parent","").split(",")[0] or i
    tx=defaultdict(lambda: dict(iv=[],keys=[],gene=None,seqid=None,strand=None))
    for f in lines:
        a,fmt=attrs(f[8])
        if fmt=="gff3" and f[2] in ("mRNA","transcript"):
            t=a.get("ID","")
            if t:
                for v in (a.get("proteinId",""),a.get("protein_id",""),a.get("Name",""),t):
                    for k in cands(v):
                        if k not in tx[t]["keys"]: tx[t]["keys"].append(k)
    for f in lines:
        if f[2]!="CDS": continue
        a,fmt=attrs(f[8])
        if fmt=="gff3":
            t=a.get("Parent","").split(",")[0]; g=m2g.get(t,t)
            pv=[a.get("protein_id",""),a.get("proteinId",""),a.get("ID",""),t]
        else:
            t=a.get("transcript_id",""); g=a.get("gene_id",t); pv=[t]
        if not t: continue
        d=tx[t]; d["gene"]=g; d["seqid"]=f[0]; d["strand"]=f[6]
        d["iv"].append((int(f[3]),int(f[4])))
        for v in pv:
            for k in cands(v):
                if k not in d["keys"]: d["keys"].append(k)
    best={}
    for t,d in tx.items():
        L=sum(e-s+1 for s,e in d["iv"])
        g=d["gene"]
        if g not in best or L>best[g][0]:
            best[g]=(L,dict(seqid=d["seqid"],strand=d["strand"],
                            iv=sorted(d["iv"]),cdslen=L,nex=len(d["iv"]),keys=d["keys"]))
    return {g:v[1] for g,v in best.items()}

def read_fa(p):
    d={}; n=None; b=[]
    with op(p) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if n: d[n]="".join(b)
                n=line[1:].split()[0]; b=[]
            else: b.append(line.strip())
    if n: d[n]="".join(b)
    return {k:v.upper().replace("*","").replace("-","").replace(".","").strip()
            for k,v in d.items()}

def attach(genes, fa):
    from collections import defaultdict as _dd
    tmp=_dd(set)
    for pid in fa:
        for k in cands(pid): tmp[k].add(pid)
    idx={k:next(iter(v)) for k,v in tmp.items() if len(v)==1}
    for g,d in genes.items():
        d["pid"]=next((idx[k] for k in d["keys"] if k in idx), None)
        d["seq"]=fa.get(d["pid"],"") if d["pid"] else ""
    return genes

def ovl_bp(a,b):
    i=j=t=0
    while i<len(a) and j<len(b):
        s=max(a[i][0],b[j][0]); e=min(a[i][1],b[j][1])
        if e>=s: t+=e-s+1
        if a[i][1]<b[j][1]: i+=1
        else: j+=1
    return t

def find(p,x):
    while p[x]!=x: p[x]=p[p[x]]; x=p[x]
    return x

def med(x): return round(st.median(x),1) if x else ""
def pc(a,b): return round(100.0*a/b,2) if b else 0.0

rows=[]
for rf in sorted(REF_FA.glob("*_longest_isoform.faa")):
    s=rf.name[:-len("_longest_isoform.faa")]
    of=OWN_FA/f"{s}_longest_isoform.faa"
    gf=next((GFF/f"{s}{e}" for e in (".gff3",".gff") if (GFF/f"{s}{e}").is_file()),None)
    tf=next((GTF/f"{s}{e}" for e in (".gtf",)  if (GTF/f"{s}{e}").is_file()),None)
    if not (of.is_file() and gf and tf): print(f"[缺文件] {s}"); continue

    R=attach(load_ann(gf), read_fa(rf))
    O=attach(load_ann(tf), read_fa(of))
    sr={d["seqid"] for d in R.values()}; so={d["seqid"] for d in O.values()}
    share=len(sr&so)
    gR=sum(1 for d in R.values() if d["seqid"] in so)
    gO=sum(1 for d in O.values() if d["seqid"] in sr)
    if share==0 or gR < 0.5*len(R):
        print(f"{s:30s} !! 序列名不匹配 (共有 {share} 个) 无法坐标配对；"
              f"参考示例 {sorted(sr)[:2]}  本方示例 {sorted(so)[:2]}")
        continue

    buck=defaultdict(list)
    for g,d in R.items(): buck[(d["seqid"],d["strand"])].append(("R",g,d))
    for g,d in O.items(): buck[(d["seqid"],d["strand"])].append(("O",g,d))
    edges=[]
    for key,lst in buck.items():
        lst.sort(key=lambda x:x[2]["iv"][0][0])
        for i,(ta,ga,da) in enumerate(lst):
            ea=da["iv"][-1][1]
            for tb,gb,db in lst[i+1:]:
                if db["iv"][0][0]>ea: break
                if ta==tb: continue
                v=ovl_bp(da["iv"],db["iv"])
                if v>=MIN_OVL_BP and v>=MIN_OVL_FRAC*min(da["cdslen"],db["cdslen"]):
                    edges.append((("R",ga) if ta=="R" else ("O",ga),
                                  ("R",gb) if tb=="R" else ("O",gb)))
    nodes=[("R",g) for g in R]+[("O",g) for g in O]
    par={n:n for n in nodes}
    for a,b in edges:
        ra,rb=find(par,a),find(par,b)
        if ra!=rb: par[ra]=rb
    comp=defaultdict(lambda:[[],[]])
    for n in nodes:
        c=find(par,n); comp[c][0 if n[0]=="R" else 1].append(n[1])

    one=split=merge=cplx=ronly=oonly=0
    exact=0; lrat=[]; ident=[]
    ro_len=[];ro_ex=[];oo_len=[];oo_ex=[]
    for c,(rs,os_) in comp.items():
        if len(rs)==1 and len(os_)==1:
            one+=1
            a,b=R[rs[0]],O[os_[0]]
            if a["seq"] and b["seq"]:
                if hashlib.md5(a["seq"].encode()).digest()==hashlib.md5(b["seq"].encode()).digest():
                    exact+=1
                lrat.append(len(b["seq"])/max(1,len(a["seq"])))
                if a["seq"] and b["seq"]:
                    K=5
                    ka={a["seq"][i:i+K] for i in range(len(a["seq"])-K+1)}
                    kb={b["seq"][i:i+K] for i in range(len(b["seq"])-K+1)}
                    if ka and kb:
                        ident.append(100.0*len(ka&kb)/len(ka|kb))
        elif len(rs)==1 and len(os_)>1: split+=1
        elif len(rs)>1 and len(os_)==1: merge+=1
        elif len(rs)>1 and len(os_)>1: cplx+=1
        elif len(rs)==1 and not os_:
            ronly+=1; ro_len.append(len(R[rs[0]]["seq"])); ro_ex.append(R[rs[0]]["nex"])
        elif len(os_)==1 and not rs:
            oonly+=1; oo_len.append(len(O[os_[0]]["seq"])); oo_ex.append(O[os_[0]]["nex"])
    allR=[len(d["seq"]) for d in R.values()]; allRe=[d["nex"] for d in R.values()]

    r=dict(sample=s, n_ref=len(R), n_own=len(O), shared_seqids=share,
        one2one=one, one2one_pct_ref=pc(one,len(R)),
        exact=exact, exact_pct_ref=pc(exact,len(R)),
        identity_median=med(ident), lenratio_median=round(st.median(lrat),4) if lrat else "",
        split=split, merge=merge, complex=cplx,
        ref_only=ronly, ref_only_pct=pc(ronly,len(R)),
        own_only=oonly, own_only_pct=pc(oonly,len(O)),
        refonly_medlen=med(ro_len), refonly_mono_pct=pc(sum(1 for x in ro_ex if x==1),len(ro_ex)),
        ref_all_medlen=med(allR),  ref_all_mono_pct=pc(sum(1 for x in allRe if x==1),len(allRe)),
        ownonly_medlen=med(oo_len), ownonly_mono_pct=pc(sum(1 for x in oo_ex if x==1),len(oo_ex)))
    rows.append(r)
    print(f"{s:30s} 1:1 {one:6d}({r['one2one_pct_ref']:5.1f}%) 精确{r['exact_pct_ref']:5.1f}% "
          f"拆{split:4d} 并{merge:4d} 杂{cplx:4d} | ref独有{ronly:6d} "
          f"中位{r['refonly_medlen']}aa/单外显子{r['refonly_mono_pct']}% "
          f"(全体 {r['ref_all_medlen']}aa/{r['ref_all_mono_pct']}%) | own独有{oonly:5d}")

if rows:
    cols=list(rows[0].keys())
    with (OUT/"coord_comparison_summary.tsv").open("w") as fh:
        fh.write("\t".join(cols)+"\n")
        for r in rows: fh.write("\t".join(str(r[c]) for c in cols)+"\n")
    print(f"\n{len(rows)} 株 -> {OUT/'coord_comparison_summary.tsv'}")
