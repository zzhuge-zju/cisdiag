
# -*- coding: utf-8 -*-
"""距离稀疏化敏感性（修正版）"""
import os
import sys, os, json
from pathlib import Path
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
CIS=os.environ.get("CISDIAG_CIS") or str(Path(__file__).resolve().parents[1] / "cis_regions")
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")
plan=json.load(open(os.path.join(D,"c2v2","_dl_plan.json")))["sel"]
cis=json.load(open(os.path.join(D,"c2v2","_cis_final.json"),encoding="utf-8"))
full={(x["a"],x["b"]):x for x in json.load(open(os.path.join(D,"c2v2","_full_region.json"),encoding="utf-8"))}
diff=[r for r in cis if plan.get(r["a"],"?").split("/")[-1]!=plan.get(r["b"],"?").split("/")[-1]]
def load(g):
    p=os.path.join(CIS,g+".tsv"); out={}
    for line in open(p,encoding="utf-8"):
        a=line.rstrip("\n").split("\t")
        if len(a)<12: continue
        try: out[(a[0].replace("chr",""),int(float(a[1])))]=(a[4],a[5],float(a[6]),float(a[7]))
        except Exception: pass
    return out
def cos_thin(A,B,th,spacing):
    rows=[]
    for k in set(A)&set(B):
        ea,oa,ba,pa=A[k]; eb,ob,bb,pb=B[k]
        if ea==eb and oa==ob: pass
        elif ea==ob and oa==eb: bb=-bb
        else: continue
        if not (pa<th or pb<th): continue      # 修正：只保留显著
        rows.append((k[1],ba,bb))
    if len(rows)<3: return None,len(rows)
    rows.sort()
    keep=[]; last=-10**12
    for pos,x,y in rows:
        if spacing>0 and pos-last<spacing: continue
        keep.append((x,y)); last=pos
    if len(keep)<3: return None,len(keep)
    M=np.array(keep); na=np.linalg.norm(M[:,0]); nb=np.linalg.norm(M[:,1])
    if na==0 or nb==0: return None,len(keep)
    return float(M[:,0]@M[:,1]/(na*nb)), len(keep)
print(f"{'pair':22s} {'full':>7s} {'10kb':>7s} {'50kb':>7s} {'100kb':>7s}  n_full n_100")
res={}
for r in sorted(diff,key=lambda x:-abs(full[(x["a"],x["b"])]["cos"])):
    A=load(r["a"]); B=load(r["b"])
    vals=[]
    for sp in [0,10000,50000,100000]:
        c,n=cos_thin(A,B,5e-8,sp); vals.append((c,n))
    lab=r["a"][:9]+"/"+r["b"][:9]
    form=" ".join((f"{abs(c):.3f}".rjust(7) if c is not None else "      -") for c,n in vals)
    print(f"{lab:22s} {form} {vals[0][1]:>6d} {vals[3][1]:>5d}")
    res[lab]={"full":vals[0][0],"s10":vals[1][0],"s50":vals[2][0],"s100":vals[3][0]}
json.dump(res,open(os.path.join(D,"c2v2","_ld_thin.json"),"w",encoding="utf-8"),ensure_ascii=False,indent=1)
for key in ["full","s10","s50","s100"]:
    k=sum(1 for v in res.values() if v[key] is not None and abs(v[key])>=0.99)
    print(f"|cos|>=0.99 在 {key}: {k}/20")
