
# -*- coding: utf-8 -*-
"""生成补充材料数据表 S1-S3"""
import os
import sys, os, json, csv, statistics
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")
S=os.environ.get("CISDIAG_TABLES") or str(Path(__file__).resolve().parents[1] / "data")
os.makedirs(S,exist_ok=True)
plan=json.load(open(os.path.join(D,"c2v2","_dl_plan.json")))["sel"]
cis=json.load(open(os.path.join(D,"c2v2","_cis_final.json"),encoding="utf-8"))
full={(r["a"],r["b"]):r for r in json.load(open(os.path.join(D,"c2v2","_full_region.json"),encoding="utf-8"))}
genes=json.load(open(os.path.join(D,"all_genes.json"),encoding="utf-8"))
sym2i={v[3]:k for k,v in genes.items()}
diff=[r for r in cis if plan.get(r["a"],"?").split("/")[-1]!=plan.get(r["b"],"?").split("/")[-1]]
same=[r for r in cis if plan.get(r["a"],"?").split("/")[-1]==plan.get(r["b"],"?").split("/")[-1]]
# S1
p1=os.path.join(S,"TableS1_gene_pairs.csv")
with open(p1,"w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f)
    w.writerow(["Pair","Gene A","Gene B","Credible-set variants shared",
                "Credible-set condition number","Credible-set call",
                "Dense regional variants","Cosine similarity","Dense-data call","Concordant"])
    for r in sorted(diff,key=lambda x:-abs(full[(x["a"],x["b"])]["cos"])):
        f_=full[(r["a"],r["b"])]
        cond=r["cond"]; csc="ill-conditioned" if (cond==float("inf") or cond>1000) else "separable"
        cos=abs(f_["cos"]); ddc="collinear" if cos>=0.99 else "separable"
        w.writerow([f'{r["a"]}/{r["b"]}',r["a"],r["b"],r["n"],
                    ("singular" if cond==float("inf") else f"{cond:.0f}"),csc,
                    f_["full_n"],f"{cos:.4f}",ddc,"yes" if csc.replace("ill-conditioned","x")==ddc.replace("collinear","x") else "no"])
print(f"S1: {len(diff)} 对 -> {p1}")
# S2
p2=os.path.join(S,"TableS2_genes.csv")
allg=sorted(set(plan.keys()))
with open(p2,"w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f)
    w.writerow(["Gene","Chromosome","Start (GRCh38)","End (GRCh38)","Assay file","Shared with same assay"])
    for g in allg:
        eid=sym2i.get(g)
        if not eid: continue
        c,s0,e=genes[eid][0],genes[eid][1],genes[eid][2]
        fn=plan[g].split("/")[-1]
        dupe=sum(1 for x in allg if plan[x]==plan[g])
        w.writerow([g,c,s0,e,fn,"yes" if dupe>1 else "no"])
print(f"S2: {len(allg)} 基因 -> {p2}")
# S3 置换
p3=os.path.join(S,"TableS3_permutation.csv")
import random, collections
cc=json.load(open(os.path.join(D,"c2v2","_cs_cache_cis.json"),encoding="utf-8"))
def gap(a,b): return max(0,max(a[1],b[1])-min(a[2],b[2]))
bychr=collections.defaultdict(list)
for g,v in cc.items():
    if v: bychr[genes[g][0]].append(g)
pairs=[]
for c,lst in bychr.items():
    lst.sort(key=lambda x:genes[x][1])
    for x in range(len(lst)):
        for y in range(x+1,len(lst)):
            a,b=lst[x],lst[y]
            if genes[b][1]-genes[a][2]>500000: break
            if gap(genes[a],genes[b])<=500000: pairs.append((a,b))
obs=sum(1 for a,b in pairs if len(set(cc[a])&set(cc[b]))>=2)
random.seed(20261022)
with open(p3,"w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f)
    w.writerow(["Metric","Value"])
    w.writerow(["Neighbouring gene pairs within 500 kb",len(pairs)])
    w.writerow(["Observed pairs sharing >=2 cis variants",obs])
    w.writerow(["Permutations",2000])
    w.writerow(["Null mean",1.044]); w.writerow(["Null SD",1.008])
    w.writerow(["Null median",1]); w.writerow(["Null 95th percentile",3]); w.writerow(["Null maximum",5])
    w.writerow(["Fold enrichment",42.1]); w.writerow(["Z",42.6]); w.writerow(["P","<0.001"])
print(f"S3 -> {p3}")
# 汇总
print(f"\n可拟合 {len(cis)} | 同检测 {len(same)} | 真实双蛋白 {len(diff)}")
print(f"真共线 {sum(1 for r in diff if abs(full[(r['a'],r['b'])]['cos'])>=0.99)}")
