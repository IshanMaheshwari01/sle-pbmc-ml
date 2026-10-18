# Lupus from a blood sample: a confounding-aware single-cell analysis

**Ishan Maheshwari** · [github.com/IshanMaheshwari01/sle-pbmc-ml](https://github.com/IshanMaheshwari01/sle-pbmc-ml)

**Question.** Can immune cells in a blood sample distinguish systemic lupus erythematosus (SLE) from health, and does that still hold once technical and demographic confounders cannot explain it?

**Data.** 1,263,676 single PBMCs from 162 SLE patients and 99 healthy controls (Perez et al., *Science* 2022), multiplexed 10x scRNA-seq over 4 processing batches. An audit found disease nearly aliased with batch (batch 1: 0 SLE / 28 healthy; batch 2: 88 / 1), skewed by ancestry, and with 15 of 17 men having SLE. Analysis cohort: 239 women (143 SLE, 96 healthy).

**Approach**
- Memory-bounded streaming of a 12 GB `.h5ad` into pseudobulk counts per sample × batch × cell type over **all** cells, plus a 300-cell-per-sample subsample for single-cell work.
- Single cell: per-batch QC; Harmony vs scVI integration, benchmarked with scib-metrics **and** a custom check that integration did not erase disease signal; CellTypist annotation check.
- Statistics (R): propeller-style composition test; pseudobulk limma-voom with `duplicateCorrelation` for technical replicates, DESeq2 as a second method; Hallmark pathways with cameraPR.
- Machine learning: L1 logistic regression, random forest and gradient boosting on composition + pseudobulk expression. Nested 5×5 cross-validation grouped by donor, with gene selection inside folds, compared against batch-only, composition-only and interferon-only baselines.

**Key results**
- **Pooled AUROC 0.99, but batch + ancestry alone give 0.87,** and the features identify a healthy donor's batch perfectly. Pooled accuracy overstates the biology.
- **Within the balanced batch (52 SLE vs 44 healthy), AUROC 0.92 (95% CI 0.86–0.97)**; shuffled labels give 0.49 (permutation p = 0.0099).
- <!-- POOL_SENTENCE_START -->**Same batch, same 10x library:** comparing only SLE and healthy samples run in the same library, AUROC is **0.97** (95% CI 0.95-0.98), so the signal holds.<!-- POOL_SENTENCE_END -->
- **Type I/II interferon response is the top pathway in all nine cell types** (FDR < 10⁻⁴), consistent across limma and DESeq2 (fold-change ρ 0.86–0.97).
- **A third of patients are interferon-low.** An interferon score cannot find them (AUROC 0.55); the full model can (0.86), using composition shifts (CD4 T cells ↓, pDCs ↓, classical monocytes ↑) and non-interferon genes, such as MAIT-cell markers down in CD8 T cells.

**Limitations.** Single cohort, with no external validation. Treatment is unrecorded, so some signal may reflect medication. The model detects established disease, not early disease.

**Tools.** Python (scanpy, scvi-tools/PyTorch on GPU, harmonypy, scib-metrics, CellTypist, scikit-learn), R/Bioconductor (limma, edgeR, DESeq2), HDF5 streaming, Make, pytest with planted-effect mock data for every pipeline step.
