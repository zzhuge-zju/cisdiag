
# -*- coding: utf-8 -*-
"""2000 次置换，给出稳健零分布"""
import os
import sys, os, json, random, collections, statistics
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')
D=os.environ.get("CISDIAG_DATA") or str(Path(__file__).resolve().parents[1] / "data")
cc=json.load(open(os.path.join(D,"c2v2","_cs_cache_cis.json"),encoding="utf-8"))
genes=json.load(open(os.path.join(D,"all_genes.json"),encoding="utf-8"))
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
obs=sum(1 for a,b in pairs if len(set(cc[a])&set(cc[b]))>=2)
random.seed(20261022)
null=[]
for it in range(2000):
    perm=[]
    for c,lst in bychr.items():
        sh=lst[:]; random.shuffle(sh)
        for i in range(len(sh)):
            for j in range(i+1,len(sh)): perm.append((sh[i],sh[j]))
    s=random.sample(perm,len(pairs))
    null.append(sum(1 for a,b in s if len(set(cc[a])&set(cc[b]))>=2))
m=statistics.mean(null); sd=statistics.stdev(null)
qs=sorted(null)
print(f"2000 次置换")
print(f"  零分布: 均值 {m:.3f} | SD {sd:.3f} | 中位 {statistics.median(null):.0f} | 95%分位 {qs[int(0.95*len(qs))]} | 最大 {max(null)}")
print(f"  观测 {obs}")
print(f"  富集 {obs/m:.1f} 倍 | Z = {(obs-m)/sd:.1f}")
print(f"  P = {sum(1 for x in null if x>=obs)/len(null):.5f}")
