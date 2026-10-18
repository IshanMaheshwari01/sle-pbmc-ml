# Phase 3: donor-level lupus classifier

- Cohort: 239 donors (143 SLE), 301 sample x batch units; same cohort as phase 2
- Features: 11 cell-type proportions (CLR); 91,539 expressed genes across 7 cell types, of which the 1000 most variable per cell type are chosen inside each training fold
- CV: 5-fold, grouped by donor, stratified, 5 repeats; L1 penalty tuned in an inner 3-fold grouped CV (chosen C: {0.3: 11, 1.0: 8, 0.1: 6})
- AUROC: 0.5 = chance, 1.0 = perfect. Brackets: 95% bootstrap CI over donors.

## Headline: accuracy pooled vs where batch and ancestry cannot help

| model                              | All donors (pooled)   | Batch 2.0 only   | Batch 3.0 only   | Batch 4.0 only   | Asian donors only   | European American donors only   |
|:-----------------------------------|:----------------------|:-----------------|:-----------------|:-----------------|:--------------------|:--------------------------------|
| Covariates only (batch + ancestry) | 0.87 [0.82-0.91]      | n/a              | n/a              | n/a              | 0.81 [0.72-0.89]    | 0.84 [0.77-0.91]                |
| Cell-type proportions only         | 0.90 [0.86-0.94]      | 0.98 [0.96-1.00] | 0.91 [0.80-0.98] | 0.84 [0.76-0.92] | 0.93 [0.87-0.97]    | 0.87 [0.81-0.93]                |
| Interferon score only              | 0.88 [0.84-0.92]      | 0.83 [0.76-0.90] | 0.95 [0.83-1.00] | 0.78 [0.68-0.86] | 0.91 [0.85-0.96]    | 0.84 [0.77-0.90]                |
| L1 logistic (all features)         | 0.99 [0.97-0.99]      | 1.00 [1.00-1.00] | 1.00 [0.97-1.00] | 0.92 [0.86-0.97] | 0.98 [0.96-1.00]    | 0.99 [0.97-1.00]                |
| Random forest (all features)       | 0.98 [0.96-0.99]      | 0.99 [0.98-1.00] | 0.99 [0.97-1.00] | 0.89 [0.82-0.95] | 0.97 [0.94-0.99]    | 0.98 [0.96-0.99]                |
| Gradient boosting (all features)   | 0.99 [0.98-1.00]      | 1.00 [0.99-1.00] | 1.00 [1.00-1.00] | 0.93 [0.87-0.98] | 0.99 [0.96-1.00]    | 0.99 [0.97-1.00]                |

`n/a`: the covariate-only model is not scored within a batch, where batch is constant; refitting a model with no usable signal fold by fold pushes its pooled held-out AUROC below 0.5, an artefact rather than information.

Read across a row: the pooled number can lean on batch (batch 1 is all healthy, batch 2 mostly lupus). Within one batch, batch is constant and cannot separate anyone, and within one ancestry group the same holds for ancestry. A model that keeps its accuracy in those columns has learned biology; the covariate-only row shows what batch and ancestry give for free.

## Stability across CV repeats (pooled, donor level)

|                                    |   mean AUROC |   SD over repeats |
|:-----------------------------------|-------------:|------------------:|
| Covariates only (batch + ancestry) |        0.88  |             0.003 |
| Cell-type proportions only         |        0.902 |             0.003 |
| Interferon score only              |        0.883 |             0.003 |
| L1 logistic (all features)         |        0.978 |             0.005 |
| Random forest (all features)       |        0.976 |             0.002 |
| Gradient boosting (all features)   |        0.984 |             0.003 |

## How much batch is written into the features?

Among **healthy donors only**, the main model's features separate batch 1.0 (46 units) from batch 4.0 (44 units) with AUROC **1.00** (same grouped CV). Batch is strongly encoded in expression, so pooled accuracy alone would overstate what the model knows about lupus.

## Permutation test (main model)

Disease labels were shuffled across donors 100 times and the whole CV re-run (penalty fixed at C = 0.3).
- Pooled: observed AUROC 0.99; null mean 0.50 (95th percentile 0.58); permutation p = 0.0099
- Batch 4.0 only: observed 0.92; null mean 0.49 (95th percentile 0.60); p = 0.0099
- With 100 permutations the smallest possible p is 0.0099.

## Lupus patients without a high interferon signature

An SLE donor counts as IFN-high when the classical-monocyte IFN score is more than 2.0 SD above healthy donors of the same batch: **96 IFN-high, 47 IFN-low** SLE donors. The table shows how well each model separates each group from healthy donors.

| model                              |   Batch 4.0 only: SLE IFN-high |   Batch 4.0 only: SLE IFN-low |   pooled: SLE IFN-high |   pooled: SLE IFN-low |
|:-----------------------------------|-------------------------------:|------------------------------:|-----------------------:|----------------------:|
| Covariates only (batch + ancestry) |                           0.37 |                          0.41 |                   0.81 |                  0.79 |
| Cell-type proportions only         |                           0.9  |                          0.78 |                   0.95 |                  0.85 |
| Interferon score only              |                           0.99 |                          0.55 |                   0.99 |                  0.63 |
| L1 logistic (all features)         |                           0.98 |                          0.86 |                   1    |                  0.97 |
| Random forest (all features)       |                           1    |                          0.76 |                   1    |                  0.92 |
| Gradient boosting (all features)   |                           1    |                          0.85 |                   1    |                  0.96 |

The IFN-only baseline should be near-perfect for IFN-high patients by construction. The IFN-low columns are the real test: anything above chance there comes from biology other than the interferon signature.

## What the main model relies on

Held-out AUROC drop when one block of features is shuffled (larger = more relied on):

|                                            |   mean |   std |
|:-------------------------------------------|-------:|------:|
| classical monocyte expression              |  0.071 | 0.03  |
| CD4-positive, alpha-beta T cell expression |  0.031 | 0.023 |
| B cell expression                          |  0.012 | 0.008 |
| CD8-positive, alpha-beta T cell expression |  0.009 | 0.011 |
| cell-type proportions                      |  0.004 | 0.004 |
| natural killer cell expression             |  0.001 | 0.004 |
| conventional dendritic cell expression     |  0.001 | 0.008 |
| non-classical monocyte expression          | -0     | 0.003 |

Features most consistently selected by the L1 model across folds (of 25 fits; positive coefficient = higher in SLE):

|                                         |   selected in % of folds |   mean coefficient |
|:----------------------------------------|-------------------------:|-------------------:|
| proportion: classical monocyte          |                      100 |              0.29  |
| B cell|EPSTI1                           |                      100 |              0.286 |
| CD4-positive, alpha-beta T cell|EPSTI1  |                      100 |              0.204 |
| classical monocyte|CXCL10               |                      100 |              0.134 |
| classical monocyte|PTCH2                |                      100 |             -0.236 |
| classical monocyte|SGK1                 |                      100 |             -0.332 |
| CD4-positive, alpha-beta T cell|RNF175  |                      100 |             -0.362 |
| conventional dendritic cell|SNAI1       |                       96 |              0.24  |
| classical monocyte|ID1                  |                       96 |              0.19  |
| classical monocyte|TNFRSF12A            |                       96 |             -0.207 |
| CD8-positive, alpha-beta T cell|EPSTI1  |                       92 |              0.196 |
| classical monocyte|ENSG00000259001      |                       92 |             -0.159 |
| CD4-positive, alpha-beta T cell|CD70    |                       88 |              0.219 |
| B cell|KLK1                             |                       84 |             -0.148 |
| CD8-positive, alpha-beta T cell|SLC4A10 |                       84 |             -0.212 |
| classical monocyte|RMI1                 |                       80 |              0.208 |
| CD8-positive, alpha-beta T cell|MCOLN2  |                       80 |              0.115 |
| classical monocyte|IL2RB                |                       80 |             -0.122 |
| classical monocyte|BLK                  |                       76 |             -0.156 |
| classical monocyte|CA2                  |                       72 |              0.221 |
| classical monocyte|IGF2BP3              |                       72 |              0.143 |
| CD4-positive, alpha-beta T cell|GBP1    |                       72 |              0.11  |
| CD4-positive, alpha-beta T cell|NUAK2   |                       72 |             -0.073 |
| B cell|SH2D1A                           |                       72 |             -0.263 |
| CD8-positive, alpha-beta T cell|CCR6    |                       68 |             -0.138 |

## Figures

- `09_auroc_by_setting.png`
- `09_roc_within_batch.png`
- `09_feature_block_importance.png`
- `09_predictions_by_ifn_status.png`

Total time: 14.7 min
