# Phase 1: cell-type label check with CellTypist

- Model: `Immune_All_Low.pkl`; majority voting over `leiden_scvi` clusters (13 clusters)
- Cells: 78,898; distinct CellTypist labels: 11
- **Overall lineage agreement with the author labels: 88.2%**

| author label                    |   cells |   lineage agreement % | top CellTypist labels                                                                              |
|:--------------------------------|--------:|----------------------:|:---------------------------------------------------------------------------------------------------|
| CD4-positive, alpha-beta T cell |   23642 |                  81.7 | Tcm/Naive helper T cells (82%); Tcm/Naive cytotoxic T cells (17%); Tem/Trm cytotoxic T cells (1%)  |
| classical monocyte              |   18891 |                  91.5 | Classical monocytes (91%); Non-classical monocytes (8%); DC2 (0%)                                  |
| CD8-positive, alpha-beta T cell |   15993 |                  85.2 | Tem/Trm cytotoxic T cells (55%); Tcm/Naive cytotoxic T cells (31%); Tcm/Naive helper T cells (13%) |
| B cell                          |    9486 |                  99.9 | B cells (100%); Tcm/Naive helper T cells (0%); Regulatory T cells (0%)                             |
| natural killer cell             |    5824 |                  96.9 | CD16+ NK cells (97%); Tem/Trm cytotoxic T cells (3%); Tcm/Naive helper T cells (0%)                |
| non-classical monocyte          |    2947 |                  98.9 | Non-classical monocytes (99%); Classical monocytes (1%); CD16+ NK cells (0%)                       |
| conventional dendritic cell     |    1088 |                  85.9 | DC2 (86%); Classical monocytes (13%); Non-classical monocytes (1%)                                 |
| lymphocyte                      |     546 |                   0   | Regulatory T cells (66%); Tem/Trm cytotoxic T cells (21%); Tcm/Naive helper T cells (7%)           |
| plasmacytoid dendritic cell     |     350 |                 100   | pDC (100%)                                                                                         |
| plasmablast                     |      77 |                   0   | Regulatory T cells (99%); B cells (1%)                                                             |
| progenitor cell                 |      54 |                   0   | Megakaryocyte precursor (100%)                                                                     |

Low agreement for a label means either the author label is coarse or ambiguous (e.g. `lymphocyte`, which has no single lineage), or the two methods genuinely disagree; the top-labels column shows which. Fine-grained CellTypist subsets (e.g. regulatory vs effector-memory helper T cells) can be used for the cell-type proportion analysis in phase 2.

Figure: `05_celltypist_vs_author.png`; per-cell predictions: `05_celltypist_predictions.csv`
