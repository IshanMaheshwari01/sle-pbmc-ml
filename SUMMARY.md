# Lupus from a blood sample: one-page summary

Ishan Maheshwari · [github.com/IshanMaheshwari01/sle-pbmc-ml](https://github.com/IshanMaheshwari01/sle-pbmc-ml)

**Question.** Can the immune cells in a blood sample tell systemic lupus erythematosus (SLE) apart from health? And does that still work once batch, ancestry and other confounders can't explain it?

**Data.** 1,263,676 single PBMCs from 162 SLE patients and 99 healthy controls (Perez et al., Science 2022), multiplexed 10x scRNA-seq over 4 processing batches. When I checked the design, disease was almost the same as batch (batch 1: 0 SLE / 28 healthy, batch 2: 88 / 1), ancestry was unbalanced, and 15 of the 17 men had SLE. So I used 239 women (143 SLE, 96 healthy).

**What I did**
- Read the 12 GB `.h5ad` file in blocks so it fits in laptop memory, and made pseudobulk counts per sample, batch and cell type from all cells, plus a 300-cell-per-sample subsample for the single-cell steps.
- Single cell: QC per batch, Harmony vs scVI integration scored with scib-metrics, plus my own check that integration didn't remove the disease signal, and a CellTypist annotation check.
- Statistics in R: propeller-style composition test, pseudobulk limma-voom with `duplicateCorrelation` for technical replicates, DESeq2 as a second method, and Hallmark pathways with cameraPR.
- Machine learning: L1 logistic regression, random forest and gradient boosting on composition and pseudobulk expression. Nested 5x5 cross-validation grouped by donor, with gene selection inside the folds, compared against batch-only, composition-only and interferon-only baselines.

**Main results**
- Pooled AUROC is 0.99, but batch and ancestry alone give 0.87, and the features can identify a healthy donor's batch perfectly. So the pooled number overstates the biology.
- In the balanced batch (52 SLE vs 44 healthy), AUROC is 0.92 (95% CI 0.86 to 0.97). Shuffled labels give 0.49 (permutation p = 0.0099).
- <!-- POOL_SENTENCE_START -->Same batch, same 10x library: comparing only SLE and healthy samples run in the same library, AUROC is 0.97 (95% CI 0.95-0.98), so the signal holds.<!-- POOL_SENTENCE_END -->
- Type I and II interferon response is the top pathway in all nine cell types (FDR < 1e-4), and limma and DESeq2 agree (fold-change correlation 0.86 to 0.97).
- A third of patients have low interferon. An interferon score can't find them (AUROC 0.55) but the full model can (0.86), using changes in cell mix (fewer CD4 T cells and pDCs, more classical monocytes) and genes outside the interferon pathway, like MAIT cell markers going down in CD8 T cells.

**Limitations.** One cohort with no external validation. Treatment isn't recorded, so some of the signal may come from medication. The model detects established disease, not early disease.

**Tools.** Python (scanpy, scvi-tools and PyTorch on GPU, harmonypy, scib-metrics, CellTypist, scikit-learn), R and Bioconductor (limma, edgeR, DESeq2), HDF5, Make, and pytest with fake data that has known effects built in.
