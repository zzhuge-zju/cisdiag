
# -*- coding: utf-8 -*-
"""cis-only 完整统计"""
import os
import sys, os, json, collections, statistics
from pathlib import Path
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")
cc=json.load(open(os.path.join(D,"c2v2","_cs_cache_cis.json"),encoding="utf-8"))
genes=json.load(open(os.path.join(D,"all_genes.json"),encoding="utf-8"))
t2d=json.load(open(os.path.join(D,"t2d_hits_v2.json"),encoding="utf-8"))
sym=lambda g: genes.get(g,("","","","?"))[3]
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
rows=[]
for a,b in sh2:
    X=[];y=[];w=[];seA=[]
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
        X.append([va["beta"],vb["beta"]]); y.append(bo); w.append(1/seo**2); seA.append(va.get("se"))
    n=len(X)
    if n<2: continue
    X=np.array(X);y=np.array(y);w=np.array(w);W=np.diag(w);M=X.T@W@X
    ev=np.linalg.eigvalsh(M)
    cond=float(ev.max()/ev.min()) if ev.min()>0 else float("inf")
    rho=None; fA=None; att=None
    if cond<1e12:
        inv=np.linalg.inv(M); coef=inv@(X.T@W@y)
        uni=float((X[:,0]@W@y)/(X[:,0]@W@X[:,0]))
        if uni!=0: att=float((coef[0]-uni)/uni*100)
        if n>=3 and np.std(X[:,0])>0:
            rho=float(np.corrcoef(X[:,0],X[:,1])[0,1])
            if all(s is not None for s in seA):
                A=np.column_stack([np.ones(n),X[:,1]])
                c=np.linalg.solve(A.T@W@A,A.T@W@X[:,0]); r_=X[:,0]-A@c
                fA=float(np.mean((r_/np.array(seA))**2))
    rows.append({"a":sym(a),"b":sym(b),"n":n,"cond":cond,"rho":rho,"fA":fA,"att":att})
json.dump(rows,open(os.path.join(D,"c2v2","_cis_final.json"),"w",encoding="utf-8"),ensure_ascii=False,indent=1)
print(f"配对 {len(rows)}")
print(f"\n{'配对':24s} {'n':>2s} {'rho':>7s} {'F':>8s} {'条件数':>12s}")
for r in sorted(rows,key=lambda x:-x["n"]):
    print(f"{r['a'][:10]+'/'+r['b'][:10]:24s} {r['n']:>2d} "
          f"{(f'{r['rho']:+.2f}' if r['rho'] is not None else '   -'):>7s} "
          f"{(f'{r['fA']:.1f}' if r['fA'] is not None else '   -'):>8s} {r['cond']:>12.0f}")
sing=sum(1 for r in rows if r["cond"]==float("inf") or r["cond"]>1e12)
big=sum(1 for r in rows if 1000<r["cond"]<1e12)
print(f"\n奇异 {sing} | 条件数>1000(非奇异) {big} | 合计病态 {sing+big}/{len(rows)}")
