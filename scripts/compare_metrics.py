
# -*- coding: utf-8 -*-
"""区分有效检验(不同检测) vs 无效检验(同一检测)，量化可信集的双向误判"""
import os
import sys, os, json
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")
plan=json.load(open(os.path.join(D,"c2v2","_dl_plan.json")))["sel"]
cis=json.load(open(os.path.join(D,"c2v2","_cis_final.json"),encoding="utf-8"))
full={ (r["a"],r["b"]):r for r in json.load(open(os.path.join(D,"c2v2","_full_region.json"),encoding="utf-8")) }
same=[];diff=[]
for r in cis:
    ka=plan.get(r["a"],"?").split("/")[-1]; kb=plan.get(r["b"],"?").split("/")[-1]
    (same if ka==kb else diff).append(r)
print(f"可拟合 29 对 → 同检测 {len(same)} 对 | 不同检测 {len(diff)} 对\n")
print("=== 不同检测配对：可信集口径 vs 全量口径 ===")
print(f"{'配对':22s} {'cs_n':>4s} {'cs_cond':>10s} {'cs判定':>7s} {'全量|cos|':>9s} {'全量判定':>7s}")
fp=0; tp=0; tn=0
for r in diff:
    f=full.get((r["a"],r["b"]),{})
    csc=r["cond"]; csbad = (csc==float("inf") or csc>1000)
    fc=abs(f.get("cos",0)) if "cos" in f else None
    fbad = fc is not None and fc>=0.99
    if csbad and not fbad: fp+=1
    if csbad and fbad: tp+=1
    if not csbad and not fbad: tn+=1
    print(f"{r['a'][:9]+'/'+r['b'][:9]:22s} {r['n']:>4d} {csc:>10.0f} {'病态' if csbad else '可分':>7s} "
          f"{(f'{fc:.4f}' if fc is not None else '   -'):>9s} {'共线' if fbad else '可分':>7s}")
print(f"\n可信集判为病态、全量实为可分（假阳性）: {fp}")
print(f"两个口径都判为不可分（真阳性）: {tp}")
print(f"两个口径都判为可分（真阴性）: {tn}")
