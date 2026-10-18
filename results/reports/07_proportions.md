# Phase 2: cell-type proportions

- Units: 301 sample x batch units from 239 donors (143 SLE, 96 healthy); cell counts from all cells
- Model: limma on logit(proportion) ~ disease_groupSLE + batch2 + batch3 + batch4 + ancestryEuropean American
- Within-donor correlation between technical replicates: 0.85

| cell type | mean % healthy | mean % SLE | odds ratio | 95% CI low | 95% CI high | FDR |
|---|---|---|---|---|---|---|
| CD4-positive, alpha-beta T cell | 37.6000 | 25.5000 | 0.542 | 0.463 | 0.635 | < 1e-4 |
| classical monocyte | 17.8000 | 28.5000 | 1.800 | 1.530 | 2.130 | < 1e-4 |
| lymphocyte | 0.3900 | 0.8430 | 1.810 | 1.480 | 2.220 | < 1e-4 |
| plasmacytoid dendritic cell | 0.6100 | 0.3120 | 0.562 | 0.451 | 0.701 | < 1e-4 |
| non-classical monocyte | 3.1100 | 4.2000 | 1.360 | 1.130 | 1.630 | 0.0021 |
| B cell | 12.6000 | 12.3000 | 0.680 | 0.525 | 0.880 | 0.0065 |
| conventional dendritic cell | 1.5800 | 1.3700 | 0.773 | 0.644 | 0.928 | 0.0092 |
| natural killer cell | 7.8800 | 6.6600 | 0.766 | 0.607 | 0.967 | 0.035 |
| progenitor cell | 0.0672 | 0.0586 | 0.731 | 0.543 | 0.984 | 0.048 |
| plasmablast | 0.0999 | 0.1330 | 0.724 | 0.518 | 1.010 | 0.063 |
| CD8-positive, alpha-beta T cell | 18.2000 | 20.1000 | 1.010 | 0.874 | 1.170 | 0.87 |

Odds ratio > 1 means that cell type makes up a larger share of a lupus donor's PBMCs.
Proportions are compositional: if one cell type expands, others shrink as a share even when
their absolute numbers do not change, so results are read together, not one by one.

Figure: `07_proportions_forest.png`
