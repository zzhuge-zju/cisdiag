
# -*- coding: utf-8 -*-
"""CI（Wilson 手算）+ Fisher + 阈值敏感性"""
import os
import sys, os, json, math
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')
from scipy.stats import fisher_exact, binomtest
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")
def wilson(k,n,z=1.96):
    if n==0: return (0,0)
    p=k/n; d=1+z*z/n
    c=p+z*z/(2*n); m=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))
    return max(0,(c-m)/d), min(1,(c+m)/d)
plan=json.load(open(os.path.join(D,"c2v2","_dl_plan.json")))["sel"]
cis=json.load(open(os.path.join(D,"c2v2","_cis_final.json"),encoding="utf-8"))
full={(r["a"],r["b"]):r for r in json.load(open(os.path.join(D,"c2v2","_full_region.json"),encoding="utf-8"))}
diff=[r for r in cis if plan.get(r["a"],"?").split("/")[-1]!=plan.get(r["b"],"?").split("/")[-1]]
vals=[abs(full[(r["a"],r["b"])]["cos"]) for r in diff]
out={}
print("=== 阈值敏感性 ===")
for th in [0.999,0.99,0.95,0.90,0.80]:
    k=sum(1 for v in vals if v>=th); lo,hi=wilson(k,20)
    print(f"  |cos| >= {th:<5}: {k}/20 = {100*k/20:3.0f}%  (95%CI {100*lo:.0f}-{100*hi:.0f}%)")
    out[f"th_{th}"]={"k":k,"ci":[round(100*lo,1),round(100*hi,1)]}
print("\n=== 关键比例 95%CI ===")
for k,n,lab in [(2,20,"真共线"),(6,20,"误判"),(9,29,"同检测"),(44,1840,"共享>=2")]:
    lo,hi=wilson(k,n); print(f"  {lab:10s} {k}/{n}: {100*k/n:5.1f}% (95%CI {100*lo:.1f}-{100*hi:.1f}%)")
    out[lab]={"k":k,"n":n,"ci":[round(100*lo,1),round(100*hi,1)]}
odds,p=fisher_exact([[1,5],[1,13]])
print(f"\nFisher 1/5/1/13: OR={odds:.2f}, P={p:.3f}")
out["fisher"]={"OR":round(float(odds),2),"P":round(float(p),3)}
r2=binomtest(2,20,0.5)
print(f"2/20 vs 50%: P={r2.pvalue:.2e}")
out["binom_2of20"]=f"{r2.pvalue:.2e}"
json.dump(out,open(os.path.join(D,"c2v2","_stats_supp.json"),"w"),indent=1)
