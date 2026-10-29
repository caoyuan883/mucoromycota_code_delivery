#!/usr/bin/env python3
import sys, re, gzip, statistics as st
from pathlib import Path
from collections import defaultdict

def op(p): return gzip.open(p,"rt") if str(p).endswith(".gz") else open(p,"r")

def attrs(s):
    s=s.strip()
    if "=" in s and not re.search(r'\w+\s+"', s):          # GFF3
        d={}
        for kv in s.rstrip(";").split(";"):
            kv=kv.strip()
            if "=" in kv:
                k,v=kv.split("=",1); d[k.strip()]=v.strip()
        return d,"gff3"
    d={}                                                    # GTF
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', s):
        d[m.group(1)]=m.group(2)
    return d,"gtf"

def collect(path):
    mrna2gene={}; cds=defaultdict(list)
    lines=[]
    with op(path) as fh:
        for line in fh:
            if not line or line[0]=="#": continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9: continue
            a,fmt=attrs(f[8]); lines.append((f[2],int(f[3]),int(f[4]),a,fmt))
    for t,s,e,a,fmt in lines:
        if fmt=="gff3" and t in ("mRNA","transcript"):
            fid=a.get("ID",""); 
            if fid: mrna2gene[fid]=a.get("Parent","").split(",")[0] or fid
    for t,s,e,a,fmt in lines:
        if t!="CDS": continue
        if fmt=="gff3":
            tid=a.get("Parent","").split(",")[0]
            gid=mrna2gene.get(tid,tid)
        else:
            tid=a.get("transcript_id",""); gid=a.get("gene_id",tid)
        if not tid: continue
        cds[(gid,tid)].append((s,e))
    return cds

def stats(cds):
    per_gene=defaultdict(list)
    for (g,t),iv in cds.items():
        per_gene[g].append((sum(e-s+1 for s,e in iv), t, sorted(iv)))
    n_genes=len(per_gene); n_tx=len(cds)
    mono=multi=0; nex=[]; cdslen=[]; intr=[]; maxex=0
    for g,lst in per_gene.items():
        lst.sort(key=lambda x:(-x[0], x[1]))
        L,t,iv=lst[0]
        k=len(iv); nex.append(k); cdslen.append(L); maxex=max(maxex,k)
        if k==1: mono+=1
        else:
            multi+=1
            for i in range(1,k): intr.append(iv[i][0]-iv[i-1][1]-1)
    med=lambda x: st.median(x) if x else 0
    return dict(n_genes=n_genes, n_transcripts=n_tx,
        isoforms_per_gene=round(n_tx/n_genes,4) if n_genes else 0,
        mono=mono, multi=multi,
        mono_multi_ratio=round(mono/multi,4) if multi else "",
        median_exons=med(nex), mean_exons=round(sum(nex)/len(nex),3) if nex else 0,
        max_exons=maxex, median_cds_nt=med(cdslen),
        median_intron_nt=med(intr), n_introns=len(intr))

if len(sys.argv)<4:
    sys.exit("用法: gene_structure_stats.py <目录> <后缀> <输出tsv>")
d=Path(sys.argv[1]); sufs=[x for x in sys.argv[2].split(",") if x]; out=Path(sys.argv[3])
cols=["n_genes","n_transcripts","isoforms_per_gene","mono","multi",
      "mono_multi_ratio","median_exons","mean_exons","max_exons",
      "median_cds_nt","median_intron_nt","n_introns"]
found={}
for sf in sufs:
    for q in sorted(d.glob("*"+sf)):
        found.setdefault(q.name[:-len(sf)], q)
rows=[]
for name in sorted(found):
    p=found[name]
    r=stats(collect(p)); rows.append((name,r))
    print(f"{name:32s} 基因{r['n_genes']:6d} 转录本{r['n_transcripts']:6d} "
          f"mono:multi={r['mono_multi_ratio']} 外显子中位{r['median_exons']}")
with out.open("w") as fh:
    fh.write("sample\t"+"\t".join(cols)+"\n")
    for n,r in rows: fh.write(n+"\t"+"\t".join(str(r[c]) for c in cols)+"\n")
print(f"\n{len(rows)} 株 -> {out}")
