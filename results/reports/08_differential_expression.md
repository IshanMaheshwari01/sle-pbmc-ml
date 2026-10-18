# Phase 2: pseudobulk differential expression

- Cohort: 239 donors, profiles with >= 10 cells; one sample per donor, technical replicates across batches kept
- limma-voom ~ disease + batch + ancestry, duplicateCorrelation on donor; robust eBayes; BH FDR < 0.05
- DESeq2 check: one profile per sample, same covariates

| cell type | SLE donors | healthy donors | genes tested | DE up | DE down | up, |lfc|>=0.5 | down, |lfc|>=0.5 | limma vs DESeq2 rho | % limma DE also DESeq2 DE |
|---|---|---|---|---|---|---|---|---|---|
| CD4-positive, alpha-beta T cell | 143 | 96 | 9430 | 1520 | 1560 | 133 | 52 | 0.97 | 82.1 |
| CD8-positive, alpha-beta T cell | 143 | 96 | 8530 | 1230 | 974 | 101 | 115 | 0.95 | 70.7 |
| classical monocyte | 143 | 96 | 9310 | 1390 | 1700 | 168 | 93 | 0.95 | 73.8 |
| B cell | 141 | 96 | 7270 | 412 | 422 | 60 | 20 | 0.91 | 67.5 |
| non-classical monocyte | 141 | 95 | 5460 | 562 | 491 | 112 | 23 | 0.91 | 68.6 |
| natural killer cell | 139 | 96 | 5720 | 326 | 362 | 43 | 8 | 0.92 | 71.5 |
| conventional dendritic cell | 136 | 95 | 4310 | 336 | 335 | 53 | 6 | 0.92 | 68.3 |
| lymphocyte | 128 | 75 | 2020 | 72 | 81 | 29 | 3 | 0.86 | 68.0 |
| plasmacytoid dendritic cell | 80 | 74 | 1380 | 103 | 116 | 27 | 2 | 0.95 | 75.8 |

## Strongest genes per cell type

| cell_type | top_up | top_down |
|---|---|---|
| CD4-positive, alpha-beta T cell | PSME2, PSMB9, CIAO2B, CARD16, TAP1, GBP1 | RPL21, RPL35A, DPEP2, RPS27A, CDKN2AIP, RPS6 |
| CD8-positive, alpha-beta T cell | PSME2, PSMB8, GBP1, ISG15, STAT1, PARP9 | SLC4A10, ATF7IP2, NR1D1, ZBTB16, RAB33B-AS1, SESN1 |
| classical monocyte | NMI, CHMP5, LAP3, PLAC8, STAT1, ANXA2 | KIF22, KRT10, EIF4B, PABPC4, PABPC1, DDX3X |
| B cell | SAMD9L, EIF2AK2, PSME2, DRAP1, STAT1, EPSTI1 | KLK1, IL10RA, NAP1L1, COBLL1, MT-CYB, BCL2 |
| non-classical monocyte | MT2A, CHMP5, GLRX, PSMB9, TNFSF10, SP110 | RPS8, KRT10, CIRBP, RPL6, NAP1L1, RPS6 |
| natural killer cell | ISG15, PARP9, PSMB9, BST2, STAT1, TAP1 | CIRBP, MGAT4A, EIF1, RPS27A, ATF7IP2, RPS14 |
| conventional dendritic cell | PSME2, LAP3, MT2A, TYMP, PSMB9, STAT1 | RPS6, RPS8, RPL6, RPL5, RPS14, RPS27A |
| lymphocyte | LY6E, BST2, IFI6, XAF1, ISG15, NKG7 | LEF1, RPS27A, RPS4X, RPL3, RPL13A, RPL7 |
| plasmacytoid dendritic cell | PSMB9, JUN, PLAC8, HSPA8, DYNLL1, ACTB | RPS27A, CXXC5, RPL32, CIRBP, RPL41, RPL35 |

## Interferon genes (log2 fold change, SLE vs healthy)

The six genes of the phase-1 interferon score, now tested properly at the donor level.

| cell_type | ISG15 | IFI27 | IFI44L | RSAD2 | IFIT1 | SIGLEC1 |
|---|---|---|---|---|---|---|
| CD4-positive, alpha-beta T cell | 1.10 | 2.20 | 1.6 | 1.2 | 1.10 | NA |
| CD8-positive, alpha-beta T cell | 0.91 | 2.30 | 1.7 | 1.4 | 0.64 | NA |
| classical monocyte | 1.40 | 2.20 | 1.7 | 1.6 | 1.70 | 1.8 |
| B cell | 0.76 | 0.98 | 1.8 | 1.1 | 1.00 | NA |
| non-classical monocyte | 1.40 | 2.90 | 1.8 | 1.6 | 1.80 | NA |
| natural killer cell | 0.91 | 1.10 | 1.9 | 1.2 | NA | NA |
| conventional dendritic cell | 1.20 | 1.40 | 1.7 | NA | 1.20 | NA |
| lymphocyte | 0.96 | 2.70 | 1.6 | NA | NA | NA |
| plasmacytoid dendritic cell | 1.30 | NA | 1.5 | NA | NA | NA |

## Pathways (h.all.symbols.gmt), FDR < 0.05

| cell_type | pathway | NGenes | Direction | FDR |
|---|---|---|---|---|
| classical monocyte | HALLMARK_INTERFERON_ALPHA_RESPONSE | 89 | Up | < 1e-4 |
| non-classical monocyte | HALLMARK_INTERFERON_ALPHA_RESPONSE | 76 | Up | < 1e-4 |
| B cell | HALLMARK_INTERFERON_ALPHA_RESPONSE | 77 | Up | < 1e-4 |
| natural killer cell | HALLMARK_INTERFERON_ALPHA_RESPONSE | 72 | Up | < 1e-4 |
| conventional dendritic cell | HALLMARK_INTERFERON_ALPHA_RESPONSE | 61 | Up | < 1e-4 |
| non-classical monocyte | HALLMARK_INTERFERON_GAMMA_RESPONSE | 143 | Up | < 1e-4 |
| classical monocyte | HALLMARK_INTERFERON_GAMMA_RESPONSE | 175 | Up | < 1e-4 |
| CD4-positive, alpha-beta T cell | HALLMARK_INTERFERON_ALPHA_RESPONSE | 85 | Up | < 1e-4 |
| CD8-positive, alpha-beta T cell | HALLMARK_INTERFERON_ALPHA_RESPONSE | 84 | Up | < 1e-4 |
| lymphocyte | HALLMARK_INTERFERON_ALPHA_RESPONSE | 43 | Up | < 1e-4 |
| CD4-positive, alpha-beta T cell | HALLMARK_INTERFERON_GAMMA_RESPONSE | 163 | Up | < 1e-4 |
| natural killer cell | HALLMARK_INTERFERON_GAMMA_RESPONSE | 133 | Up | < 1e-4 |
| lymphocyte | HALLMARK_INTERFERON_GAMMA_RESPONSE | 67 | Up | < 1e-4 |
| CD8-positive, alpha-beta T cell | HALLMARK_INTERFERON_GAMMA_RESPONSE | 163 | Up | < 1e-4 |
| B cell | HALLMARK_INTERFERON_GAMMA_RESPONSE | 151 | Up | < 1e-4 |
| conventional dendritic cell | HALLMARK_INTERFERON_GAMMA_RESPONSE | 120 | Up | < 1e-4 |
| plasmacytoid dendritic cell | HALLMARK_INTERFERON_ALPHA_RESPONSE | 28 | Up | < 1e-4 |
| plasmacytoid dendritic cell | HALLMARK_INTERFERON_GAMMA_RESPONSE | 48 | Up | < 1e-4 |
| non-classical monocyte | HALLMARK_COMPLEMENT | 105 | Up | < 1e-4 |
| CD4-positive, alpha-beta T cell | HALLMARK_IL2_STAT5_SIGNALING | 130 | Up | 0.00011 |
| CD8-positive, alpha-beta T cell | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 194 | Up | 0.00018 |
| CD4-positive, alpha-beta T cell | HALLMARK_INFLAMMATORY_RESPONSE | 100 | Up | 0.00044 |
| plasmacytoid dendritic cell | HALLMARK_HYPOXIA | 27 | Up | 0.00047 |
| classical monocyte | HALLMARK_COMPLEMENT | 136 | Up | 0.001 |
| classical monocyte | HALLMARK_TNFA_SIGNALING_VIA_NFKB | 163 | Down | 0.001 |

## Files

- Per-gene tables: `08_de/<cell type>_limma.csv` (and `_deseq2.csv`)
- Figures: `08_de_counts.png`, `08_volcano_classical_monocyte.png`, `08_limma_vs_deseq2.png`, `08_pathways_heatmap.png`
