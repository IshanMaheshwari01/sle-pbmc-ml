# Can a blood sample tell you who has lupus?

I used 1.26 million single blood cells from 162 lupus (SLE) patients and 99 healthy people to see whether immune cells can separate the two groups. The catch is that in this dataset, disease status is badly mixed up with processing batch and ancestry. So the real question became: how much of the accuracy is left once batch, ancestry and 10x library can't explain it?

## Short answer

- If you just pool everyone, the classifier gets an AUROC of 0.99. But batch and ancestry alone already get 0.87, and the same features can tell which batch a healthy donor was processed in with AUROC 1.00. So the pooled number is inflated.
- Inside batch 4, the only batch with similar numbers of lupus and healthy donors (52 vs 44), batch can't help. There the model gets an AUROC of 0.92 (95% CI 0.86 to 0.97). With shuffled labels it gets 0.49 (permutation p = 0.0099, the smallest possible with 100 permutations). I treat 0.92 as the real result.
- It isn't explained by the 10x library either. Comparing only lupus and healthy samples that were run in the same library gives AUROC 0.97 overall and 0.92 in batch 4.
- The interferon signature shows up in every immune cell type. Interferon alpha and gamma response are the top pathways in all nine cell types I tested (FDR < 1e-4), after adjusting for batch and ancestry, and limma-voom and DESeq2 agree (fold-change correlation 0.86 to 0.97).
- The model does more than detect interferon. About a third of the patients (47 of 143) have low interferon. In batch 4, an interferon score alone can't separate them from healthy donors (AUROC 0.55), but the full model can (0.86).
- Lupus blood has fewer CD4 T cells (odds ratio 0.54) and plasmacytoid DCs (0.56), and more classical monocytes (1.80). In the CD8 T cell group, the genes that drop the most are MAIT cell markers (SLC4A10, ZBTB16, KLRB1, CCR6), which suggests MAIT cells are lost.

![Classifier AUROC pooled and within batch / ancestry](results/figures/09_auroc_by_setting.png)

## The problem with this dataset

The data is from Perez et al. 2022 (Science): 1,263,676 PBMCs, multiplexed 10x 3' v2, in 23 pools across 4 processing batches. Before doing any analysis I checked how the groups were spread out (`results/reports/01_data_audit.md`), and disease status was tangled up with almost everything:

| Variable | What I found | What I did about it |
|---|---|---|
| Processing batch | Batch 1: 0 lupus / 28 healthy. Batch 2: 88 / 1. 54 samples were sequenced in two batches | Batch is in every model. I evaluate within batches that have both groups, keep cross-batch replicates as separate profiles, and compare against a batch-only baseline |
| Ancestry | Asian donors: 83 lupus / 24 healthy. European: 75 / 74 | Ancestry is in every model, I evaluate per ancestry, and I left out 5 donors from very small groups |
| Sex | 15 of the 17 men have lupus | Women only (239 donors: 143 lupus, 96 healthy) |
| 10x pool | Several donors share each library | Within-library evaluation (section 5) |
| Repeated samples | Technical replicates, plus flare and post-treatment samples from the same patient | One sample per donor, all CV splits grouped by donor, and `duplicateCorrelation` for replicates |

## Results

### 1. Interferon signalling is up across the immune system

![Hallmark pathways by cell type](results/figures/08_pathways_heatmap.png)

I ran pseudobulk limma-voom for each cell type with `expression ~ disease + batch + ancestry`, and used `duplicateCorrelation` on donor for replicates. Interferon alpha and gamma response are at the top in every cell type. The six interferon score genes (IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1) are 2 to 7 times higher in lupus, in T cells, B cells, NK cells and monocytes. DESeq2 on one profile per sample agrees: fold-change correlation 0.86 to 0.97, and 68 to 82% of the genes significant in limma are also significant in DESeq2. Tables are in `results/reports/08_differential_expression.md`.

### 2. Lupus blood has a different mix of cells

![Cell type proportions](results/figures/07_proportions_forest.png)

For proportions I used all 1.26M cells and modelled each sample on the logit scale with limma (the propeller approach). The drop in CD4 T cells and pDCs and the rise in monocytes match the original paper. The B cell result (odds ratio 0.68) only appears after adjusting for batch, because the raw percentages are almost the same, so I'd call that one model-dependent.

### 3. The classifier, and what's left once batch can't help

| Model | Pooled | Batch 4 only (52 vs 44) | Batch 3 only (15 vs 16) | Batch 2 only (106 vs 22) |
|---|---|---|---|---|
| Batch + ancestry only | 0.87 | n/a | n/a | n/a |
| Interferon score only | 0.88 | 0.78 | 0.95 | 0.83 |
| Cell type proportions only | 0.90 | 0.84 | 0.91 | 0.98 |
| L1 logistic, all features | 0.99 | 0.92 | 1.00 | 1.00 |
| Gradient boosting, all features | 0.99 | 0.93 | 1.00 | 1.00 |
| Random forest, all features | 0.98 | 0.89 | 0.99 | 0.99 |

- Features: cell type composition (centred log-ratio) and pseudobulk log-CPM for seven cell types. The 1,000 most variable genes per cell type are picked inside each training fold, so the test fold never influences gene selection.
- Validation: 5-fold CV grouped by donor, repeated 5 times. The L1 penalty is tuned with an inner grouped CV.

Batch 4 is the one I trust most because it's the only balanced batch. The near-perfect scores in batches 2 and 3 aren't caused by 10x pools either (section 5). My guess for why batch 4 is harder: its lupus patients have the weakest interferon signature of the three batches (median monocyte IFN score 0.03, vs 0.12 in batch 2 and 0.46 in batch 3), which could mean milder or better-controlled disease.

### 4. Finding lupus without the interferon signature

![Predictions by interferon status](results/figures/09_predictions_by_ifn_status.png)

I called a lupus donor interferon-high if their monocyte interferon score was more than 2 SD above the healthy donors in the same batch. That gives 96 IFN-high and 47 IFN-low donors. In batch 4, the interferon score can't tell IFN-low patients from healthy donors (0.55, basically chance), but the full model gets 0.86. A lot of that comes from cell composition (0.78 from proportions alone), plus expression changes outside the interferon pathway, like lower MAIT cell markers in CD8 T cells and CD70 in CD4 T cells.

### 5. Pool check: same batch, same 10x library

<!-- POOL_CHECK_START -->

Comparing only lupus and healthy samples run in the same 10x library (so neither batch nor library can help), the classifier's AUROC is 0.97 (95% CI 0.95-0.98, 3,324 pairs). The signal holds: the model separates lupus from healthy samples run in the same 10x library.

| Setting | Libraries (with both groups) | Samples with a same-library partner | Cramér's V, library × diagnosis | Within-library AUROC, model | Within-library, IFN only | Ordinary AUROC, model |
|---|---|---|---|---|---|---|
| Batch 1 | 12 (0) | 0% |  |  |  |  |
| Batch 2 | 40 (36) | 93% | 0.18 | 1.00 | 0.87 | 1.00 |
| Batch 3 | 12 (12) | 100% | 0.13 | 1.00 | 0.98 | 0.98 |
| Batch 4 | 24 (24) | 100% | 0.09 | 0.92 | 0.83 | 0.92 |
| All batches | 88 (72) | 82% | 0.53 | 0.97 | 0.86 | 0.98 |

This section is filled in by `make poolcheck`. Details are in `results/reports/10_pool_check.md`.

<!-- POOL_CHECK_END -->

Within each batch the libraries are well mixed (Cramér's V 0.09 to 0.19), so most samples have a partner from the other group in the same library.

### 6. Integration and annotation

- I compared Harmony and scVI with scib-metrics. Harmony scored highest (0.72 vs 0.65 for scVI and 0.61 uncorrected), mostly because it mixed batches better.
- I also checked, within batch 4, that neither method removed the disease signal (classical monocytes: 1.25 / 1.25 / 1.29 enrichment of same-diagnosis neighbours from other donors).
- I only used the integrated embeddings for plots and annotation. All statistics and models use uncorrected pseudobulk counts with batch as a covariate, because when batch and disease are this mixed up, batch correction can remove real disease signal.
- CellTypist agreed with the authors' labels for 92% of cells. Where they disagreed it was interesting: the authors' "lymphocyte" label is mostly cycling and regulatory T cells, and a third of the conventional DCs sit with monocytes, which is the known DC2 / monocyte overlap.

## What this doesn't show

- Treatment. Most lupus donors have "managed" disease and are on medication, but the dataset doesn't record which. Some of the signal could be drug effects. For example, steroids suppress NF-kB signalling and change monocyte and lymphocyte counts.
- External validation. Everything is cross-validated inside one cohort. Testing on a second lupus cohort would be the real test.
- The ancestry comparisons still contain batch, so they're a weaker check than the within-batch results.
- Diagnosis, not prognosis. The model separates people who already have lupus from healthy people. It says nothing about early or undiagnosed disease.

## How it's built

| Step | Scripts | What it does |
|---|---|---|
| 0. Data | `00_download.py`, `01_inspect.py`, `02_subsample.py` | Downloads the 12 GB file (can resume), runs the confounding audit, then reads the file once in blocks to make a 300-cell-per-sample subsample and pseudobulk over all cells. Memory depends on block size, not file size |
| 1. Single cell | `03_qc.py`, `04_integrate.py`, `05_annotate.py` | QC per batch (MAD-based), Harmony vs scVI on the GPU with scib-metrics and the disease-signal check, interferon score, CellTypist with majority voting |
| 2. Statistics (R) | `06_export_pseudobulk.py`, `R/07_proportions.R`, `R/08_differential_expression.R` | Applies the cohort rules, propeller-style proportions, limma-voom with `duplicateCorrelation`, DESeq2 as a second check, cameraPR on MSigDB Hallmark |
| 3. Classifier | `09_classifier.py` | Nested CV grouped by donor, three baselines, AUROC within batch and per ancestry with bootstrap CIs, batch detectability, 100-run permutation test, IFN-low analysis, feature block importance |
| 4. Checks | `10_pool_check.py` | AUROC within each 10x library, and how mixed each batch's libraries are |

Every step has tests that run on a small fake dataset with the same layout as the real file. The fake data has known effects built in (an interferon signal in lupus monocytes, a change in composition, a batch effect), and the tests check that each one is found and that nothing is found where nothing was put in.

## Running it

```bash
mamba create -n sle-py python=3.11 -y && conda activate sle-py
pip install -e ".[phase1]" pytest            # scvi-tools brings PyTorch with CUDA
sudo apt install -y r-base-core r-bioc-limma r-bioc-edger r-bioc-deseq2 r-cran-ggplot2 r-cran-yaml r-cran-matrix

make test         # 26 tests on fake data, about 3 min
make download     # 12.2 GB from CELLxGENE, can resume
make inspect      # confounding audit, read this first
make genesets     # MSigDB Hallmark gene sets
make all          # subsample, step 1, step 2, classifier, pool check
```

I ran it on a laptop with 7.4 GB RAM available to WSL and an RTX 4050 (6 GB). Peak RAM is about 5 GB. Rough times: subsample 1.5 min, step 1 about 25 min (scVI takes 10 min on the GPU), step 2 about 1 hour (mostly DESeq2), classifier 15 min, pool check a few minutes.

## Folders

```
config/       config.yaml (paths, columns, subsampling), phase1/2/3.yaml (analysis settings), mock.yaml
scripts/      00 to 10, one Python script per step
R/            step 2 statistics
src/sle_pbmc/ h5ad block reader, sampling and pseudobulk, helpers for steps 1 and 3, plot style
tests/        fake data builder and tests for every step
results/      reports/ (markdown and csv for every step) and figures/
envs/         conda environments
```

## Data and citation

Perez RK, Gordon MG, Subramaniam M, et al. Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. Science 376, eabf1970 (2022). Data from [CELLxGENE](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2) (also GEO GSE174188). Genotypes are controlled-access and I didn't use them.

A one-page version of this is in [SUMMARY.md](SUMMARY.md).
