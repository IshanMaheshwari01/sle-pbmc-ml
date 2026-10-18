# Phase 1: integration

- Cells: 78,898; highly variable genes: 3,000 (seurat_v3, ranked across `Processing_Cohort`); PCs: 50
- scVI: 400 epochs on NVIDIA GeForce RTX 4050 Laptop GPU in 10.5 min; latent dimensions: 30
- Leiden clusters on the scVI graph: 16

## Standard integration benchmark (scib-metrics)

_computed on 24,998 cells, stratified by batch and cell type_

|                    | type             |   Unintegrated (PCA) |   Harmony |   scVI |
|:-------------------|:-----------------|---------------------:|----------:|-------:|
| Isolated labels    | Bio conservation |                0.604 |     0.606 |  0.583 |
| KMeans NMI         | Bio conservation |                0.687 |     0.72  |  0.641 |
| KMeans ARI         | Bio conservation |                0.447 |     0.506 |  0.434 |
| Silhouette label   | Bio conservation |                0.566 |     0.568 |  0.542 |
| cLISI              | Bio conservation |                0.996 |     0.997 |  0.996 |
| BRAS               | Batch correction |                0.898 |     0.94  |  0.912 |
| iLISI              | Batch correction |                0.425 |     0.607 |  0.473 |
| KBET               | Batch correction |                0.477 |     0.728 |  0.577 |
| Graph connectivity | Batch correction |                0.892 |     0.87  |  0.943 |
| PCR comparison     | Batch correction |                0     |     0.805 |  0.412 |
| Batch correction   | Aggregate score  |                0.538 |     0.79  |  0.664 |
| Bio conservation   | Aggregate score  |                0.66  |     0.68  |  0.639 |
| Total              | Aggregate score  |                0.611 |     0.724 |  0.649 |

All metrics are 0-1, higher is better (unscaled).

**Read the batch-correction scores with care.** They reward cells from different batches mixing together. Here batch 1 is all healthy and batch 2 almost all lupus, so part of what a high batch score rewards is lupus and healthy cells overlapping, i.e. lost biology. The next table checks that directly.

## Is disease signal still there? (batch 4.0 only)

For each cell, its nearest neighbours **from other donors** are taken in each embedding, within one cell type and within the one batch that holds both groups. *Enrichment* is how much more often those neighbours share the cell's disease status than random cells would (1.0 = none). If integration lowers enrichment compared with the unintegrated PCA, it removed disease differences along with batch.

| cell_type                       |   cells |   Unintegrated (PCA) |   Harmony |   scVI |
|:--------------------------------|--------:|---------------------:|----------:|-------:|
| CD4-positive, alpha-beta T cell |    8575 |                1.078 |     1.078 |  1.158 |
| classical monocyte              |    5467 |                1.248 |     1.249 |  1.287 |
| CD8-positive, alpha-beta T cell |    4056 |                1.111 |     1.111 |  1.16  |
| B cell                          |    2546 |                1.088 |     1.087 |  1.147 |
| natural killer cell             |    1883 |                1.111 |     1.11  |  1.169 |
| non-classical monocyte          |     774 |                1.334 |     1.334 |  1.331 |
| conventional dendritic cell     |     284 |                1.117 |     1.125 |  1.097 |

## Interferon score in classical monocytes

Six-gene type I IFN score (Rice et al. 2013) from log-normalised expression, median per sample, averaged per donor. AUROC = how well the score alone separates lupus from healthy donors (0.5 = chance). This uses no batch correction, so a high AUROC *within* a batch cannot come from batch.

| batch                   |   SLE donors |   healthy donors |   median SLE |   median healthy |   AUROC | p (Mann-Whitney)   |
|:------------------------|-------------:|-----------------:|-------------:|-----------------:|--------:|:-------------------|
| 1.0                     |            0 |               45 |      nan     |           -0.3   | nan     |                    |
| 2.0                     |          120 |               21 |        0.117 |           -0.26  |   0.799 | < 1e-4             |
| 3.0                     |           22 |               16 |        0.458 |           -0.258 |   0.895 | < 1e-4             |
| 4.0                     |           52 |               44 |        0.028 |           -0.243 |   0.785 | < 1e-4             |
| all (one row per donor) |          161 |               98 |        0.2   |           -0.268 |   0.84  | < 1e-4             |

Donors sequenced in more than one batch (technical replicates) appear in each of those batch rows, which is what puts healthy donors into batch 2. The pooled row mixes the disease effect with the batch effect; the within-batch rows are the ones to trust.

## Figures

- `04_umap_integration.png`
- `04_umap_cell_types.png`
- `04_ifn_score.png`

Total time: 14.4 min
