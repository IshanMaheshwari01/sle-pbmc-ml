# Can a blood sample tell you who has lupus?

**Donor-level machine learning and pseudobulk statistics on 1.26 million single PBMCs from 162 lupus patients and 99 healthy controls.**

Systemic lupus erythematosus (SLE) is diagnosed from a mix of symptoms and antibody tests, and it looks different from patient to patient. This project asks two questions of a large single-cell dataset:

1. **Can a model separate lupus patients from healthy donors using only their blood immune cells, and does it still work on samples from a processing batch it has never seen?**
2. **Which cell types carry the signal, and through which genes and pathways?**

The headline is not the accuracy number. It is whether that number survives the checks that most single-cell ML papers skip: grouping by donor, holding out whole batches, and testing whether batch, ancestry or sex alone can predict the label.

## Data

Perez et al. 2022, *Science* ([PMC9297655](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655)), downloaded from [CELLxGENE](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2) (also GEO [GSE174188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188)).

| | |
|---|---|
| Cells | 1,263,676 PBMCs after the authors' QC and doublet removal |
| People | 162 SLE patients, 99 healthy controls; European and Asian ancestry, mostly women |
| Design | Multiplexed 10x 3′ v2 in 23 pools over 4 processing batches, with technical replicates and some flare → post-treatment pairs |
| File | One 12.2 GB `.h5ad`; raw counts in `raw/X` |

Genotypes for this cohort are controlled-access (dbGaP) and are **not** used.

## Why this is harder than it looks

| Trap | What goes wrong | What this repo does |
|---|---|---|
| Cells treated as replicates | Hundreds of thousands of "samples" give tiny p-values for effects driven by a handful of people | All statistics are on **pseudobulk** profiles: one per sample × cell type |
| Same person on both sides of a CV split | Technical replicates and flare/treated pairs leak identity into the test set | Every split is **grouped by donor** |
| Batch, ancestry or sex tracks disease | The model learns the batch, not lupus | Audit in phase 0 (it does, see below); **within-batch** and **European-only** evaluation plus a **covariate-only baseline** in phase 3 |
| Thinning cells for speed | Fewer cells per sample means noisier pseudobulk | Pseudobulk uses **all** 1.26M cells; only the cell-level object is subsampled |

## What the data audit found

`make inspect` (full report in `results/reports/01_data_audit.md`) confirmed 1,263,676 cells, 261 donors (162 SLE, 99 healthy) and 274 samples, with integer raw counts. It also showed that disease status is badly confounded with other variables:

| Variable | Finding (donors) | Consequence |
|---|---|---|
| Processing batch | Batch 1: 0 SLE / 28 healthy. Batch 2: 88 SLE / 1 healthy. Only batches 3 and 4 contain both groups; 54 samples were run in more than one batch | Batch alone nearly predicts disease. Leave-one-batch-out is impossible (batches 1 and 2 have one class). Evaluation has to happen within batches that contain both groups, and cross-batch replicates are kept as separate profiles to estimate the batch effect |
| Ancestry | Asian: 83 SLE / 24 healthy. European: 75 / 74. Three African American and two Hispanic donors | Ancestry-linked expression could pass for disease signal. European donors form a nearly balanced test set; the five donors in small groups are left out of the models |
| Sex | 15 of 17 men have SLE | Sex-chromosome genes would predict disease in men. Models use female donors only (97 healthy, 147 SLE) |

The project's main question is therefore no longer just "how accurate is the classifier?" but "how much of its accuracy survives once batch, ancestry and sex can't explain it?"

## Roadmap

- [x] **Phase 0: data.** Resumable download, data audit (design, confounding, count check), one streaming pass that writes a cell-level subsample and full pseudobulk.
- [ ] **Phase 1: single-cell processing (Python).** Code ready (`make phase1`):
  - `03_qc.py`: per-batch adaptive QC check (median ± MADs). Doublet detection is not re-run, because it needs whole 10x lanes, not a 300-cell-per-sample subsample; the authors already removed doublets with freemuxlet and Scrublet.
  - `04_integrate.py`: HVGs ranked across batches, then three embeddings (uncorrected PCA, **Harmony**, **scVI**) compared with **scib-metrics** *and* with a disease-signal check, because with this much confounding a high batch-mixing score can mean lupus and healthy cells were blended together. Also a six-gene type I interferon score per donor, tested within batches.
  - `05_annotate.py`: **CellTypist** (`Immune_All_Low`) against the author labels, with lineage-level agreement and finer subsets.
  - Integrated embeddings are used only for visualisation, clustering and annotation; phases 2–3 work on pseudobulk raw counts with batch as a covariate.
- [ ] **Phase 2: statistics (R).** Code ready (`make phase2`):
  - `06_export_pseudobulk.py` writes the pseudobulk counts (from all 1.26M cells) for R and fixes the cohort once: women only, ancestry groups with at least 10 donors, one biological sample per donor.
  - `R/07_proportions.R`: cell-type proportions per sample, logit-transformed and modelled with **limma** (the propeller approach), adjusted for batch and ancestry.
  - `R/08_differential_expression.R`: pseudobulk DE per cell type with **limma-voom**, adjusted for batch and ancestry. A sample sequenced in two batches gives two profiles, handled with `duplicateCorrelation` blocked on donor. **DESeq2** is run as a second method. **Hallmark** pathways are tested with limma's `cameraPR`.
  - Tests plant a composition shift and an interferon signal in a mock and check both are recovered, and that nothing is called where nothing was planted.
- [ ] **Phase 3: donor-level ML.** Features from cell-type proportions and per-cell-type pseudobulk expression. Models: L1-logistic regression, random forest, XGBoost. Nested cross-validation with `StratifiedGroupKFold` on donor. AUROC and AUPRC with bootstrap CIs, plus calibration. Stress tests:
  - a **covariate-only baseline** (batch + ancestry, no expression): the score any expression model has to beat
  - performance **within batch 4** (25 SLE vs 41 healthy) and **within batch 3**, where batch cannot explain the prediction
  - performance among **European donors only** (75 vs 74), where ancestry cannot explain it
  - female donors only throughout; a label-permutation test for the null distribution

  Then **SHAP** to see which cell types and genes drive predictions.
- [ ] **Phase 4: write-up.** End-to-end `make all`, a Quarto report, and a one-page summary of findings and limitations.
- [ ] *Stretch:* flare vs managed disease, within lupus patients only.

## Quickstart

```bash
mamba env create -f envs/python.yml && conda activate sle-py
pip install -e .

make test        # 17 tests on a fake dataset with the real file's layout (seconds)
make mock        # run the whole phase-0 pipeline on that fake dataset

make download    # 12.2 GB from CELLxGENE; resumes if interrupted
make inspect     # writes results/reports/01_data_audit.md - read it before anything else
make dev         # tiny subsample (40 cells/sample) for writing phase-1 code quickly
make subsample   # full run: 500 cells/sample + pseudobulk over every cell
```

Phase 1 needs extra packages: `pip install -e ".[phase1]"`. Check the GPU is visible with `python -c "import torch; print(torch.cuda.is_available())"`. Run `make phase1-dev` first (tiny subsample, a few minutes) to catch problems, then `make phase1`. CellTypist downloads its model on first use.

Phase 2 needs R with limma, edgeR, DESeq2 and ggplot2 (see `envs/r.yml`; on Ubuntu/WSL one `apt install` line). Run `make genesets` once, then `make phase2`.

Have about 15 GB free disk for the raw file and outputs. `make subsample` needs roughly 5 GB of RAM at peak with the default settings (fine on a 16 GB laptop); lower `subsample.max_cells_per_sample` or `subsample.block_size` in `config/config.yaml` if memory is tight. Phase 1 (scVI) is much faster on an NVIDIA GPU; free Colab or Kaggle GPUs work if you don't have one.

## Outputs of phase 0

| File | Contents |
|---|---|
| `results/reports/01_data_audit.md` | Cells, donors, samples; disease × batch / ancestry / sex tables; technical replicates; cell-type mix |
| `results/reports/01_samples.csv` | One row per sample: donor, disease, cells, libraries, batches |
| `data/processed/cells_subsample.h5ad` | Raw counts for every sample, ≤ 500 cells each; `layers["counts"]` kept for scVI |
| `data/processed/cells_dev.h5ad` | Same, ≤ 40 cells per sample |
| `data/processed/pseudobulk_sample_celltype.h5ad` | Summed raw counts per sample × batch × cell type over all cells, with `profile_id`, `sample_id`, `donor`, covariates, `n_cells`, `n_libraries`, and a `low_cell_count` flag |

## Design notes

- **Sample definition.** A sample is `donor_id × disease_state`, so a flare and a post-treatment sample from the same patient are two samples of one donor. Set in `config/config.yaml` (`sample_key`).
- **Technical replicates.** Libraries from the same sample *within* a batch are summed. Replicates in *different* batches are kept as separate pseudobulk profiles (`pseudobulk.split_by: ["Processing_Cohort"]`), because with batch this confounded with disease, the same sample measured in two batches is the cleanest information there is about the batch effect. Each profile carries `profile_id` (sample × batch) and `sample_id` (sample); models group by donor. Set `split_by: []` to sum everything instead.
- **Streaming reader** (`src/sle_pbmc/h5ad_stream.py`). AnnData's backed mode does not handle `.raw` well, so the 12 GB file is read straight from HDF5 in contiguous row blocks. Memory depends on the block size, not the file size.
- **All column names live in the config.** If CELLxGENE renames a field in a later release, `make inspect` lists the closest matches and only the config changes.

## Repository layout

```
config/       config.yaml (real data), mock.yaml (fake data), phase1.yaml, phase2.yaml (parameters)
envs/         python.yml, r.yml
scripts/      00_download.py, 01_inspect.py, 02_subsample.py   (phase 0)
              03_qc.py, 04_integrate.py, 05_annotate.py          (phase 1)
              06_export_pseudobulk.py                            (phase 2, hands data to R)
src/sle_pbmc/ streaming h5ad reader, sampling, pseudobulk, phase-1 helpers, figure style
tests/        mock dataset builder + end-to-end tests for phases 0, 1 and 2
R/            phase 2 statistics: common.R, 07_proportions.R, 08_differential_expression.R
notebooks/    phase 3 analysis (to come)
data/         raw/ and processed/ (git-ignored)
results/      reports and figures
```

## Citation

Perez RK, Gordon MG, Subramaniam M, et al. Single-cell RNA-seq reveals cell type–specific molecular and genetic associations to lupus. *Science* 376, eabf1970 (2022).
