# Phase 1: integration

- Cells: 78,898; highly variable genes: 3,000 (seurat_v3, ranked across `Processing_Cohort`); PCs: 50
- scVI: 200 epochs on NVIDIA GeForce RTX 4050 Laptop GPU in 5.3 min; latent dimensions: 30
- Leiden clusters on the scVI graph: 13

## Standard integration benchmark (scib-metrics)

_computed on 24,998 cells, stratified by batch and cell type_

|                    | type             |   Unintegrated (PCA) |   Harmony |   scVI |
|:-------------------|:-----------------|---------------------:|----------:|-------:|
| Isolated labels    | Bio conservation |                0.604 |     0.606 |  0.577 |
| KMeans NMI         | Bio conservation |                0.687 |     0.72  |  0.712 |
| KMeans ARI         | Bio conservation |                0.447 |     0.506 |  0.512 |
| Silhouette label   | Bio conservation |                0.566 |     0.568 |  0.535 |
| cLISI              | Bio conservation |                0.996 |     0.997 |  0.993 |
| BRAS               | Batch correction |                0.898 |     0.94  |  0.92  |
| iLISI              | Batch correction |                0.425 |     0.607 |  0.478 |
| KBET               | Batch correction |                0.477 |     0.728 |  0.567 |
| Graph connectivity | Batch correction |                0.892 |     0.87  |  0.952 |
| PCR comparison     | Batch correction |                0     |     0.805 |  0.351 |
| Batch correction   | Aggregate score  |                0.538 |     0.79  |  0.654 |
| Bio conservation   | Aggregate score  |                0.66  |     0.68  |  0.666 |
| Total              | Aggregate score  |                0.611 |     0.724 |  0.661 |

All metrics are 0-1, higher is better (unscaled).

**Read the batch-correction scores with care.** They reward cells from different batches mixing together. Here batch 1 is all healthy and batch 2 almost all lupus, so part of what a high batch score rewards is lupus and healthy cells overlapping, i.e. lost biology. The next table checks that directly.

## Is disease signal still there? (batch 4.0 only)

For each cell, its nearest neighbours **from other donors** are taken in each embedding, within one cell type and within the one batch that holds both groups. *Enrichment* is how much more often those neighbours share the cell's disease status than random cells would (1.0 = none). If integration lowers enrichment compared with the unintegrated PCA, it removed disease differences along with batch.

| cell_type                       |   cells |   Unintegrated (PCA) |   Harmony |   scVI |
|:--------------------------------|--------:|---------------------:|----------:|-------:|
| CD4-positive, alpha-beta T cell |    8575 |                1.078 |     1.078 |  1.126 |
| classical monocyte              |    5467 |                1.248 |     1.249 |  1.243 |
| CD8-positive, alpha-beta T cell |    4056 |                1.111 |     1.111 |  1.129 |
| B cell                          |    2546 |                1.088 |     1.087 |  1.12  |
| natural killer cell             |    1883 |                1.111 |     1.11  |  1.124 |
| non-classical monocyte          |     774 |                1.334 |     1.334 |  1.304 |
| conventional dendritic cell     |     284 |                1.117 |     1.125 |  1.076 |

## Interferon score in classical monocytes

Six-gene type I IFN score (Rice et al. 2013) from log-normalised expression, median per sample, averaged per donor. AUROC = how well the score alone separates lupus from healthy donors (0.5 = chance). This uses no batch correction, so a high AUROC *within* a batch cannot come from batch.

| batch   |   SLE donors |   healthy donors |   median SLE |   median healthy |    AUROC |   p (Mann-Whitney) |
|:--------|-------------:|-----------------:|-------------:|-----------------:|---------:|-------------------:|
| 1.0     |            0 |               45 |     nan      |          -0.2997 | nan      |                nan |
| 2.0     |          120 |               21 |       0.1167 |          -0.2604 |   0.7992 |                  0 |
| 3.0     |           22 |               16 |       0.4578 |          -0.2579 |   0.8949 |                  0 |
| 4.0     |           52 |               44 |       0.0276 |          -0.2433 |   0.7854 |                  0 |
| all     |          194 |              126 |       0.1693 |          -0.2692 |   0.8346 |                  0 |

The `all` row pools batches and so mixes the disease effect with the batch effect; the within-batch rows are the honest ones.

## Figures

- `04_umap_integration.png`
- `04_umap_cell_types.png`
- `04_ifn_score.png`

Total time: 9.2 min
