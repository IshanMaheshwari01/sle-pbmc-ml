# Phase 1: integration

- Cells: 10,507; highly variable genes: 3,000 (seurat_v3, ranked across `Processing_Cohort`); PCs: 50
- scVI: 20 epochs on NVIDIA GeForce RTX 4050 Laptop GPU in 0.1 min; latent dimensions: 30
- Leiden clusters on the scVI graph: 17

## Standard integration benchmark (scib-metrics)

_computed on 10,507 cells, stratified by batch and cell type_

|                    | type             |   Unintegrated (PCA) |   Harmony |   scVI |
|:-------------------|:-----------------|---------------------:|----------:|-------:|
| Isolated labels    | Bio conservation |                0.59  |     0.591 |  0.348 |
| KMeans NMI         | Bio conservation |                0.708 |     0.718 |  0.708 |
| KMeans ARI         | Bio conservation |                0.509 |     0.518 |  0.514 |
| Silhouette label   | Bio conservation |                0.565 |     0.567 |  0.653 |
| cLISI              | Bio conservation |                0.994 |     0.994 |  0.997 |
| BRAS               | Batch correction |                0.859 |     0.909 |  0.726 |
| iLISI              | Batch correction |                0.441 |     0.621 |  0.468 |
| KBET               | Batch correction |                0.441 |     0.756 |  0.597 |
| Graph connectivity | Batch correction |                0.845 |     0.83  |  0.889 |
| PCR comparison     | Batch correction |                0     |     0.761 |  0.555 |
| Batch correction   | Aggregate score  |                0.517 |     0.775 |  0.647 |
| Bio conservation   | Aggregate score  |                0.673 |     0.678 |  0.644 |
| Total              | Aggregate score  |                0.611 |     0.717 |  0.645 |

All metrics are 0-1, higher is better (unscaled).

**Read the batch-correction scores with care.** They reward cells from different batches mixing together. Here batch 1 is all healthy and batch 2 almost all lupus, so part of what a high batch score rewards is lupus and healthy cells overlapping, i.e. lost biology. The next table checks that directly.

## Is disease signal still there? (batch 4.0 only)

For each cell, its nearest neighbours **from other donors** are taken in each embedding, within one cell type and within the one batch that holds both groups. *Enrichment* is how much more often those neighbours share the cell's disease status than random cells would (1.0 = none). If integration lowers enrichment compared with the unintegrated PCA, it removed disease differences along with batch.

| cell_type                       |   cells |   Unintegrated (PCA) |   Harmony |   scVI |
|:--------------------------------|--------:|---------------------:|----------:|-------:|
| CD4-positive, alpha-beta T cell |    1138 |                1.092 |     1.091 |  1.093 |
| classical monocyte              |     660 |                1.241 |     1.24  |  1.235 |
| CD8-positive, alpha-beta T cell |     569 |                1.064 |     1.067 |  1.026 |
| B cell                          |     352 |                1.079 |     1.08  |  1.013 |
| natural killer cell             |     263 |                1.074 |     1.071 |  1.058 |

## Interferon score in classical monocytes

Six-gene type I IFN score (Rice et al. 2013) from log-normalised expression, median per sample, averaged per donor. AUROC = how well the score alone separates lupus from healthy donors (0.5 = chance). This uses no batch correction, so a high AUROC *within* a batch cannot come from batch.

| batch   |   SLE donors |   healthy donors |   median SLE |   median healthy |    AUROC |   p (Mann-Whitney) |
|:--------|-------------:|-----------------:|-------------:|-----------------:|---------:|-------------------:|
| 1.0     |            0 |                5 |     nan      |          -0.2664 | nan      |           nan      |
| 2.0     |           52 |                0 |       0.2309 |         nan      | nan      |           nan      |
| 3.0     |           10 |                2 |       0.5968 |          -0.1286 | nan      |           nan      |
| 4.0     |           14 |                9 |       0.3693 |          -0.1912 |   0.8889 |             0.0022 |
| all     |           76 |               16 |       0.3147 |          -0.2459 |   0.884  |             0      |

The `all` row pools batches and so mixes the disease effect with the batch effect; the within-batch rows are the honest ones.

## Figures

- `04_umap_integration.png`
- `04_umap_cell_types.png`
- `04_ifn_score.png`

Total time: 1.7 min
