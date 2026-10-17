# Data audit: Perez et al. 2022 lupus PBMC

- File: `perez2022_sle_pbmc.h5ad`
- Cells: **1,263,676**; genes in count matrix: **30,172**
- Raw counts read from: `raw/X`
- First 200,000 stored values: all whole numbers (min 1.0, max 477.0, stored as float32)

## All obs columns

`library_uuid`, `assay_ontology_term_id`, `mapped_reference_annotation`, `is_primary_data`, `cell_type_ontology_term_id`, `author_cell_type`, `cell_state`, `sample_uuid`, `tissue_ontology_term_id`, `development_stage_ontology_term_id`, `disease_state`, `suspension_enriched_cell_types`, `suspension_uuid`, `suspension_type`, `donor_id`, `self_reported_ethnicity_ontology_term_id`, `disease_ontology_term_id`, `sex_ontology_term_id`, `Processing_Cohort`, `ct_cov`, `ind_cov`, `tissue_type`, `cell_type`, `assay`, `disease`, `sex`, `tissue`, `self_reported_ethnicity`, `development_stage`, `observation_joinid`

## Design

- Donors: **261**
- Samples (`donor_id x disease_state`): **274**

Samples per donor:

|   samples_per_donor |   donors |
|--------------------:|---------:|
|                   1 |      250 |
|                   2 |        9 |
|                   3 |        2 |

### Donors by disease

| disease                      |   count |
|:-----------------------------|--------:|
| systemic lupus erythematosus |     162 |
| normal                       |      99 |

### Disease x disease_state (samples)

| disease                      |   flare |   managed |   na |   treated |
|:-----------------------------|--------:|----------:|-----:|----------:|
| systemic lupus erythematosus |      19 |       146 |    0 |        10 |
| normal                       |       0 |         0 |   99 |         0 |

## Confounding check (donors)

Uneven rows here are covariates a classifier could exploit instead of disease. They must go into the DE model and be tested in validation.

### batch (`Processing_Cohort`) x disease

| Processing_Cohort   |   systemic lupus erythematosus |   normal |   All |
|:--------------------|-------------------------------:|---------:|------:|
| 1.0                 |                              0 |       28 |    28 |
| 2.0                 |                             88 |        1 |    89 |
| 3.0                 |                             16 |        5 |    21 |
| 4.0                 |                             25 |       41 |    66 |
| mixed               |                             33 |       24 |    57 |
| All                 |                            162 |       99 |   261 |

### ancestry (`self_reported_ethnicity`) x disease

| self_reported_ethnicity   |   systemic lupus erythematosus |   normal |   All |
|:--------------------------|-------------------------------:|---------:|------:|
| African American          |                              3 |        0 |     3 |
| European American         |                             75 |       74 |   149 |
| Hispanic or Latin         |                              1 |        1 |     2 |
| Asian                     |                             83 |       24 |   107 |
| All                       |                            162 |       99 |   261 |

### sex (`sex`) x disease

| sex    |   systemic lupus erythematosus |   normal |   All |
|:-------|-------------------------------:|---------:|------:|
| female |                            147 |       97 |   244 |
| male   |                             15 |        2 |    17 |
| All    |                            162 |       99 |   261 |

Samples spread over more than one batch (technical replicates): **54**. Their libraries are summed into one pseudobulk profile per sample, and CV folds are grouped by donor so replicates never sit on both sides of a split.

## Cells per sample

|       |   n_cells |
|:------|----------:|
| count |       274 |
| mean  |      4612 |
| std   |      2187 |
| min   |       456 |
| 25%   |      3232 |
| 50%   |      3930 |
| 75%   |      5241 |
| max   |     13543 |

## Cell types (% of cells within each disease group)

| cell_type                       |   systemic lupus erythematosus |   normal |
|:--------------------------------|-------------------------------:|---------:|
| B cell                          |                          11.82 |    12.28 |
| lymphocyte                      |                           0.82 |     0.38 |
| natural killer cell             |                           6.98 |     7.87 |
| CD4-positive, alpha-beta T cell |                          25.13 |    38.06 |
| CD8-positive, alpha-beta T cell |                          20.6  |    18.25 |
| plasmacytoid dendritic cell     |                           0.29 |     0.61 |
| classical monocyte              |                          28.46 |    17.72 |
| non-classical monocyte          |                           4.34 |     3.1  |
| plasmablast                     |                           0.12 |     0.1  |
| conventional dendritic cell     |                           1.38 |     1.54 |
| progenitor cell                 |                           0.06 |     0.07 |

Cell-level percentages are dominated by samples with many cells; the proportion analysis in phase 2 is done per sample instead.
