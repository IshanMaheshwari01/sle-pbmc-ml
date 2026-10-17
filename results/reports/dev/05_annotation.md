# Phase 1: cell-type label check with CellTypist

- Model: `Immune_All_Low.pkl`; majority voting over `leiden_scvi` clusters (17 clusters)
- Cells: 10,507; distinct CellTypist labels: 11
- **Overall lineage agreement with the author labels: 82.4%**

| author label                    |   cells |   lineage agreement % | top CellTypist labels                                                           |
|:--------------------------------|--------:|----------------------:|:--------------------------------------------------------------------------------|
| CD4-positive, alpha-beta T cell |    3118 |                  98.9 | Tcm/Naive helper T cells (99%); NK cells (1%); Tem/Trm cytotoxic T cells (0%)   |
| classical monocyte              |    2471 |                  96.9 | Classical monocytes (97%); Non-classical monocytes (3%); DC2 (0%)               |
| CD8-positive, alpha-beta T cell |    2107 |                  24.3 | Tcm/Naive helper T cells (39%); NK cells (35%); Tem/Trm cytotoxic T cells (24%) |
| B cell                          |    1306 |                  99.6 | B cells (56%); Naive B cells (43%); Tcm/Naive helper T cells (0%)               |
| natural killer cell             |     771 |                  99.5 | CD16+ NK cells (92%); NK cells (8%); Tem/Trm cytotoxic T cells (1%)             |
| non-classical monocyte          |     442 |                  97.7 | Non-classical monocytes (98%); Classical monocytes (2%)                         |
| conventional dendritic cell     |     146 |                  86.3 | DC2 (86%); Classical monocytes (14%)                                            |
| lymphocyte                      |      81 |                   0   | Cycling T cells (79%); CD16+ NK cells (9%); Tem/Trm cytotoxic T cells (6%)      |
| plasmacytoid dendritic cell     |      49 |                  95.9 | pDC (96%); Naive B cells (2%); DC2 (2%)                                         |
| progenitor cell                 |       9 |                   0   | pDC (78%); Cycling T cells (22%)                                                |
| plasmablast                     |       7 |                   0   | pDC (71%); Cycling T cells (14%); B cells (14%)                                 |

Low agreement for a label means either the author label is coarse or ambiguous (e.g. `lymphocyte`, which has no single lineage), or the two methods genuinely disagree; the top-labels column shows which. Fine-grained CellTypist subsets (e.g. regulatory vs effector-memory helper T cells) can be used for the cell-type proportion analysis in phase 2.

Figure: `05_celltypist_vs_author.png`; per-cell predictions: `05_celltypist_predictions.csv`
