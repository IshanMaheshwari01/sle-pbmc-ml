# Phase 2 cohort

- Profiles exported: 3,535 (sample x batch x cell type), genes: 30,172
- Rules: sex in ['female']; ancestry groups with >= 10 donors (kept: European American, Asian); one sample per donor, priority ['na', 'managed', 'flare', 'treated']
- Ancestry groups left out: African American (3), Hispanic or Latin (2)
- **Cohort: 239 donors (143 SLE, 96 healthy), 301 sample x batch units, 3,174 profiles**

## Donors by ancestry

| ancestry          |   SLE |   healthy |   All |
|:------------------|------:|----------:|------:|
| Asian             |    77 |        24 |   101 |
| European American |    66 |        72 |   138 |
| All               |   143 |        96 |   239 |

## SLE samples by disease state

| state   |   donors |
|:--------|---------:|
| managed |      132 |
| flare   |       11 |

## Sample x batch units by batch

| batch   |   SLE |   healthy |   All |
|:--------|------:|----------:|------:|
| 1.0     |     0 |        46 |    46 |
| 2.0     |   106 |        22 |   128 |
| 3.0     |    15 |        16 |    31 |
| 4.0     |    52 |        44 |    96 |
| All     |   173 |       128 |   301 |

A sample measured in two batches gives two units; the models treat them as correlated measurements of the same donor (limma duplicateCorrelation).
