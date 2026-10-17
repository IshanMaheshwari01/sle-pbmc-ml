# Phase 1: cell-type label check with CellTypist

- Model: `Immune_All_Low.pkl`; majority voting over 397 fine Leiden clusters on the scVI graph
- Cells: 78,898; distinct CellTypist labels: 21
- **Lineage agreement with the author labels: 91.9% after majority voting, 87.7% per cell**

| author label                    |   cells |   agreement % (majority) |   agreement % (per cell) | top CellTypist labels                                                                                 |
|:--------------------------------|--------:|-------------------------:|-------------------------:|:------------------------------------------------------------------------------------------------------|
| CD4-positive, alpha-beta T cell |   23642 |                     93.8 |                     93.5 | Tcm/Naive helper T cells (85%); Tem/Effector helper T cells (5%); Tcm/Naive cytotoxic T cells (4%)    |
| classical monocyte              |   18891 |                     98.5 |                     88.7 | Classical monocytes (99%); Non-classical monocytes (1%); DC2 (0%)                                     |
| CD8-positive, alpha-beta T cell |   15993 |                     80.3 |                     70.8 | Tcm/Naive cytotoxic T cells (29%); Tem/Temra cytotoxic T cells (26%); Tem/Trm cytotoxic T cells (19%) |
| B cell                          |    9486 |                     99.8 |                     96.4 | Naive B cells (38%); B cells (33%); Memory B cells (25%)                                              |
| natural killer cell             |    5824 |                     95.7 |                     98.4 | CD16+ NK cells (84%); NK cells (12%); Tem/Temra cytotoxic T cells (3%)                                |
| non-classical monocyte          |    2947 |                     91.7 |                     97.7 | Non-classical monocytes (92%); Classical monocytes (8%); CD16+ NK cells (0%)                          |
| conventional dendritic cell     |    1088 |                     66.6 |                     72   | DC2 (61%); Classical monocytes (33%); DC1 (5%)                                                        |
| lymphocyte                      |     546 |                      0   |                      0   | Cycling T cells (50%); Regulatory T cells (42%); Tem/Trm cytotoxic T cells (5%)                       |
| plasmacytoid dendritic cell     |     350 |                    100   |                     99.7 | pDC (100%)                                                                                            |
| plasmablast                     |      77 |                     96.1 |                     93.5 | Plasma cells (96%); Cycling T cells (3%); Memory B cells (1%)                                         |
| progenitor cell                 |      54 |                      0   |                     29.6 | Megakaryocyte precursor (98%); pDC (2%)                                                               |

Low agreement for a label means either the author label is coarse or ambiguous (e.g. `lymphocyte`, which has no single lineage), or the two methods genuinely disagree; the top-labels column shows which. Fine-grained CellTypist subsets (e.g. regulatory vs effector-memory helper T cells) can be used for the cell-type proportion analysis in phase 2.

Figure: `05_celltypist_vs_author.png`; per-cell predictions: `05_celltypist_predictions.csv`
