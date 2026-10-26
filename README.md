# cisdiag

Analysis code for **"Credible sets misjudge instrument collinearity in
cis-Mendelian randomization"** (Zhao & Zhuge).

Reproduces every number in the manuscript from the four derived tables in
`data/`.

## Quick start

The headline result needs no downloads:

```
pip install -r requirements.txt
python scripts/reproduce_headline.py
```

Output:

```
Pairs with distinct assays      : 20
Truly collinear (dense data)    : 2  (10%)
  95% CI                        : 3-30%
Flagged ill-conditioned by set  : 6
Discordant classifications      : 6  (30%)
```

## What is in `data/`

| File | Reproduces |
|---|---|
| `TableS1_gene_pairs.csv` | Figure 4 and the 2x2 comparison: the 20 gene pairs, both criteria, and whether they agree |
| `TableS2_genes.csv` | The 52 genes, their coordinates, and their source assays |
| `TableS3_permutation.csv` | The permutation control (2,000 permutations) |
| `TableS4_somascan_replication.csv` | The SomaScan replication of all 15 comparable pairs |

## Full pipeline

The scripts below require the source summary statistics, which are not
redistributed here. Set the input locations with environment variables.

| Variable | Points to | Default |
|---|---|---|
| `CISDIAG_DATA` | directory holding the derived caches and gene coordinates | `data/` |
| `CISDIAG_CIS` | directory of extracted cis-region TSVs | `cis_regions/` |

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `cis_only.py` | Open Targets credible sets | cis-only cache |
| 2 | `perm2000.py` | cis-only cache + gene coordinates | permutation null |
| 3 | `full_region.py` | UKB-PPP cis-region extractions | dense-data cosine similarity |
| 4 | `compare_metrics.py` | steps 1-3 | the 2x2 concordance table |
| 5 | `stats_supp2.py` | step 4 | confidence intervals, threshold sensitivity |
| 6 | `ld_thin3.py` | cis-region extractions | LD-thinning robustness |

## Three implementation notes

**cis filtering is mandatory.** Open Targets `credibleSets` returns both cis
and trans pQTL rows. Without filtering by chromosome and position relative to
the gene body, 56.2% of returned variants fall outside the cis window, which
inflates the number of "shared" pairs roughly four-fold. This was the largest
error in our own early iterations, and `cis_only.py` exists to prevent it.

**Same-assay pairs must be excluded.** The Olink panel assays AMY1A/B/C with
one antibody, and likewise for CKMT1A/B, FUT3/5, LGALS7/B, DEFA1/B, DEFB4A/B
and DEFB104A/B. A pair formed from these genes is one measurement, not two
proteins, and will show perfect agreement by construction. `TableS2_genes.csv`
marks them.

**Singularity detection must be deterministic.** `numpy.linalg.inv` on a
near-singular matrix can succeed or fail depending on the platform BLAS. We
test eigenvalues directly: a matrix is declared non-identifiable when the
smallest eigenvalue is <= 0 or the condition number exceeds 1e12.

## Data sources

Cis-pQTL credible sets: Open Targets Platform.
Cis-eQTL summary statistics: eQTLGen.
Outcome summary statistics: FinnGen release 12.
Dense regional pQTL summary statistics: UKB-PPP Olink release distributed by
deCODE Genetics (access request required).
LD reference: 1000 Genomes European, via LDlink.

## Licence

MIT. See `LICENSE`.
