
# -*- coding: utf-8 -*-
"""
用 UKB-PPP 全量 cis 区变异（而非仅 credible set）重跑 29 对可识别性检验
"""
import os
import sys, os, json, math
from pathlib import Path
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
CIS=os.environ.get("CISDIAG_CIS") or str(Path(__file__).resolve().parents[1] / "cis_regions")
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")

def load(g):
    p=os.path.join(CIS,g+".tsv")
    if not os.path.exists(p): return None
    out={}
    with open(p,encoding="utf-8") as f:
        for line in f:
            a=line.rstrip("\n").split("\t")
            if len(a)<12: continue
            try:
                pos=int(float(a[1])); b=float(a[6]); pv=float(a[7]); se=float(a[9])
            except: continue
            out[(a[0].replace("chr",""),pos)]=(a[4],a[5],b,se,pv)
    return out

def harmonise(xa, xb):
    """返回 (betaA, betaB) 或 None"""
    ea,oa,ba,sea,pa = xa
    eb,ob,bb,seb,pb = xb
    if ea==eb and oa==ob: return ba,bb
    if ea==ob and oa==eb: return ba,-bb
    return None

def analyse(a,b,thresh):
    A=load(a); B=load(b)
    if A is None or B is None: return None
    keys=set(A)&set(B)
    X=[]; sig=[]
    for k in keys:
        h=harmonise(A[k],B[k])
        if h is None: continue
        ba,bb=h
        pa,pb=A[k][4],B[k][4]
        if pa<thresh or pb<thresh:
            X.append([ba,bb]); sig.append(k)
    if len(X)<2: return {"n":len(X)}
    X=np.array(X)
    na=np.linalg.norm(X[:,0]); nb=np.linalg.norm(X[:,1])
    if na==0 or nb==0: return {"n":len(X)}
    cos=float(X[:,0]@X[:,1]/(na*nb))
    W=np.eye(len(X)); M=X.T@W@X
    ev=np.linalg.eigvalsh(M)
    cond=float(ev.max()/ev.min()) if ev.min()>0 else float("inf")
    rho=float(np.corrcoef(X[:,0],X[:,1])[0,1]) if len(X)>=3 and np.std(X[:,0])>0 and np.std(X[:,1])>0 else None
    return {"n":len(X),"cos":cos,"cond":cond,"rho":rho}

rows=json.load(open(os.path.join(D,"c2v2","_cis_final.json"),encoding="utf-8"))
print(f"重跑 {len(rows)} 对，用全量 cis 区显著变异（P<5e-8）\n")
print(f"{'配对':22s} {'cs_n':>5s} {'full_n':>7s} {'|cos|':>7s} {'cond':>10s} {'结论':>8s}")
out=[]
for r in rows:
    a,b=r["a"],r["b"]
    res=analyse(a,b,5e-8)
    if not res or res.get("n",0)<2:
        print(f"{a[:9]+'/'+b[:9]:22s} {r['n']:>5d} {'-':>7s}  (全量数据不足)")
        out.append({"a":a,"b":b,"cs_n":r["n"],"full_n":res.get("n",0) if res else 0,"status":"no_data"})
        continue
    coll = abs(res["cos"])>=0.99
    print(f"{a[:9]+'/'+b[:9]:22s} {r['n']:>5d} {res['n']:>7d} {abs(res['cos']):>7.4f} {res['cond']:>10.0f}  {'共线' if coll else '可分离':>8s}")
    out.append({"a":a,"b":b,"cs_n":r["n"],"full_n":res["n"],"cos":res["cos"],
                "cond":res["cond"],"rho":res.get("rho"),"collinear_full":coll})
json.dump(out,open(os.path.join(D,"c2v2","_full_region.json"),"w",encoding="utf-8"),ensure_ascii=False,indent=1)
n=len([x for x in out if "cos" in x])
c=len([x for x in out if x.get("collinear_full")])
print(f"\n可分析 {n} 对 | 全量口径下共线 {c} 对 ({100*c/n if n else 0:.0f}%)")
