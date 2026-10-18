# Phase 4: pool (10x library) check

- Cohort units: 301 (sample x batch), libraries: 88, of which 72 hold both lupus and healthy samples
- Classifier: L1 logistic (all features), re-run with 2 CV repeats (same features and folds as phase 3)

## Per batch

| setting     |   SLE units |   healthy units |   libraries |   libraries with both groups |   % units with an opposite-group partner |   Cramer's V (library x diagnosis) |   within-library AUROC (model) |   pairs compared | model 95% CI   |   ordinary AUROC (model) |   within-library AUROC (IFN only) |   within-library AUROC (proportions) |
|:------------|------------:|----------------:|------------:|-----------------------------:|-----------------------------------------:|-----------------------------------:|-------------------------------:|-----------------:|:---------------|-------------------------:|----------------------------------:|-------------------------------------:|
| Batch 1.0   |           0 |              46 |          12 |                            0 |                                    0     |                            nan     |                        nan     |                0 |                |                  nan     |                           nan     |                              nan     |
| Batch 2.0   |         106 |              22 |          40 |                           36 |                                   92.969 |                              0.185 |                          1     |             1552 | 1.00-1.00      |                    1     |                             0.869 |                                0.974 |
| Batch 3.0   |          15 |              16 |          12 |                           12 |                                  100     |                              0.125 |                          1     |              260 | 1.00-1.00      |                    0.983 |                             0.985 |                                0.854 |
| Batch 4.0   |          52 |              44 |          24 |                           24 |                                  100     |                              0.093 |                          0.923 |             1512 | 0.89-0.95      |                    0.92  |                             0.833 |                                0.844 |
| All batches |         173 |             128 |          88 |                           72 |                                   81.728 |                              0.526 |                          0.965 |             3324 | 0.95-0.98      |                    0.984 |                             0.862 |                                0.906 |

- **% units with an opposite-group partner**: share of samples that sit in a library with at least one sample of the other group. Only these can enter a within-library comparison.
- **Cramer's V (library x diagnosis)**: 0 = diagnosis unrelated to library (well mixed); 1 = every library holds only lupus or only healthy samples.
- **within-library AUROC**: the chance that a lupus sample scores above a healthy sample **from the same library**. Batch and library are identical within each comparison, so they cannot explain it.
- **ordinary AUROC**: the usual score over all units in that setting, for comparison.

## Verdict

Within-library AUROC over all batches: **0.97** (3,324 lupus-healthy pairs sharing a library; 95% CI over libraries 0.95-0.98). The classifier's signal holds: the model separates lupus from healthy samples run in the same 10x library.
