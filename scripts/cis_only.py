
# -*- coding: utf-8 -*-
"""cis-only 重建并重跑全部分析"""
import os
import sys, os, json, collections
from pathlib import Path
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")
cache=json.load(open(os.path.join(D,"c2v2","_cs_cache.json"),encoding="utf-8"))
genes=json.load(open(os.path.join(D,"all_genes.json"),encoding="utf-8"))
t2d=json.load(open(os.path.join(D,"t2d_hits_v2.json"),encoding="utf-8"))
# --- cis 过滤 ---
cc={}
for g,v in cache.items():
    if not v or g not in genes: continue
    gc,gs,ge=genes[g][0],genes[g][1],genes[g][2]
    lo,hi=gs-1_000_000,ge+1_000_000
    cc[g]={vid:x for vid,x in v.items() if str(x["chr"])==gc and lo<=x["pos"]<=hi}
json.dump(cc,open(os.path.join(D,"c2v2","_cs_cache_cis.json"),"w",encoding="utf-8"))
tot=sum(len(v) for v in cc.values())
print(f"cis-only 缓存: {sum(1 for v in cc.values() if v):,} 基因 | {tot:,} 变异")
sz=sorted(len(v) for v in cc.values() if v)
print(f"  每基因中位 {sz[len(sz)//2]} | 最大 {max(sz)}")
def gap(a,b): return max(0, max(a[1],b[1]) - min(a[2],b[2]))
bychr=collections.defaultdict(list)
for g,v in cc.items():
    if v: bychr[genes[g][0]].append(g)
pairs=[]
for c,lst in bychr.items():
    lst.sort(key=lambda g: genes[g][1])
    for i in range(len(lst)):
        for j in range(i+1,len(lst)):
            a,b=lst[i],lst[j]
            if genes[b][1]-genes[a][2]>500000: break
            if gap(genes[a],genes[b])<=500000: pairs.append((a,b))
sh2=[(a,b) for a,b in pairs if len(set(cc[a])&set(cc[b]))>=2]
print(f"\n窗口内配对: {len(pairs):,} | 共享>=2: {len(sh2)}")
ov=collections.Counter(len(set(cc[a])&set(cc[b])) for a,b in pairs)
print("重叠分布:", {k:ov[k] for k in sorted(ov)[:6]})
fit=[];bad=0;sing=0
for a,b in sh2:
    X=[];y=[];w=[]
    for vid in set(cc[a])&set(cc[b]):
        va,vb=cc[a][vid],cc[b][vid]
        if va["beta"] is None or vb["beta"] is None: continue
        o=t2d.get(f"{va['chr']}:{va['pos']}")
        if not o: continue
        if o[2]==va["ref"] and o[3]==va["alt"]: s=1
        elif o[2]==va["alt"] and o[3]==va["ref"]: s=-1
        else: continue
        try: bo=float(o[8])*s; seo=float(o[9])
        except Exception: continue
        if seo<=0: continue
        X.append([va["beta"],vb["beta"]]); y.append(bo); w.append(1/seo**2)
    n=len(X)
    if n<2: continue
    X=np.array(X);y=np.array(y);w=np.array(w);W=np.diag(w);M=X.T@W@X
    ev=np.linalg.eigvalsh(M)
    if ev.min()<=0 or ev.max()/ev.min()>1e12: sing+=1; continue
    cond=float(ev.max()/ev.min())
    if cond>1000: bad+=1
    fit.append(n)
print(f"可拟合: {len(fit)} | 奇异: {sing} | 条件数>1000: {bad}")
print(f"工具数分布: {dict(sorted(collections.Counter(fit).items()))}")
