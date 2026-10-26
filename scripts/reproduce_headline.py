#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproduce the headline numbers from the deposited tables."""
import csv, os, sys, math

def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    p = k / n; d = 1 + z*z/n
    c = p + z*z/(2*n); m = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return max(0.0, (c-m)/d), min(1.0, (c+m)/d)

def main(tables):
    rows = list(csv.DictReader(open(os.path.join(tables, "TableS1_gene_pairs.csv"), encoding="utf-8-sig")))
    n = len(rows)
    cs_bad = [r for r in rows if r["Credible-set call"] != "separable"]
    dense_bad = [r for r in rows if r["Dense-data call"] == "collinear"]
    discordant = [r for r in rows if r["Concordant"] == "no"]
    print(f"Pairs with distinct assays      : {n}")
    print(f"Truly collinear (dense data)    : {len(dense_bad)}  ({100*len(dense_bad)/n:.0f}%)")
    lo, hi = wilson(len(dense_bad), n)
    print(f"  95% CI                        : {100*lo:.0f}-{100*hi:.0f}%")
    print(f"Flagged ill-conditioned by set  : {len(cs_bad)}")
    print(f"Discordant classifications      : {len(discordant)}  ({100*len(discordant)/n:.0f}%)")
    print()
    print("Discordant pairs:")
    for r in discordant:
        print(f"  {r['Pair']:<22} credible-set={r['Credible-set call']:<16} dense |cos|={r['Cosine similarity']}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
