#!/usr/bin/env Rscript
# Phase 2, step 2: do cell-type proportions differ between lupus and healthy donors?
#
# Unit of analysis: one sample in one processing batch (cell counts come from
# all 1.26M cells, not the subsample). For each cell type the proportion is
# logit-transformed and modelled with limma, as in the propeller method
# (Phipson et al. 2022, Bioinformatics):
#     logit(proportion) ~ disease + batch + ancestry
# A sample run in two batches gives two correlated units; limma's
# duplicateCorrelation, blocked on donor, accounts for that.
#
# Usage: Rscript R/07_proportions.R [data/processed/phase2]

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "common.R"))
suppressPackageStartupMessages(library(limma))

d <- load_phase2()
p2 <- d$p2
pr <- d$profiles[d$profiles$in_cohort, ]

# units x cell types matrix of cell counts
units <- unique(pr$profile_id)
types <- sort(unique(pr$cell_type))
N <- matrix(0, nrow = length(units), ncol = length(types), dimnames = list(units, types))
N[cbind(match(pr$profile_id, units), match(pr$cell_type, types))] <- pr$n_cells
tot <- rowSums(N)
keep <- tot >= p2$proportions$min_cells_per_profile
N <- N[keep, , drop = FALSE]
tot <- tot[keep]
meta <- pr[match(rownames(N), pr$profile_id), c("profile_id", "donor", "disease_group", "batch", "ancestry")]
rownames(meta) <- meta$profile_id

prop <- N / tot
logit_p <- t(qlogis((N + 0.5) / (tot + 1)))   # cell types x units; +0.5 avoids log(0)

design <- make_design(meta, p2$model$covariates)
coef <- "disease_groupSLE"
corr <- NA
if (anyDuplicated(meta$donor)) {
  dc <- duplicateCorrelation(logit_p, design, block = meta$donor)
  corr <- dc$consensus.correlation
  fit <- lmFit(logit_p, design, block = meta$donor, correlation = corr)
} else {
  fit <- lmFit(logit_p, design)
}
fit <- eBayes(fit, robust = TRUE)
tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = TRUE)

is_sle <- meta$disease_group == "SLE"
res <- data.frame(
  cell_type = rownames(tt),
  mean_pct_healthy = 100 * colMeans(prop[!is_sle, rownames(tt), drop = FALSE]),
  mean_pct_SLE = 100 * colMeans(prop[is_sle, rownames(tt), drop = FALSE]),
  odds_ratio = exp(tt$logFC), ci_low = exp(tt$CI.L), ci_high = exp(tt$CI.R),
  t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val, check.names = FALSE
)
res <- res[order(res$p), ]
write.csv(res, file.path(d$man$reports_dir, "07_proportions.csv"), row.names = FALSE)

# ---- figure: odds ratios with 95% CI ------------------------------------------------
res$label <- short_ct(res$cell_type)
res$direction <- ifelse(res$fdr >= p2$de$fdr, "not significant",
                        ifelse(res$odds_ratio > 1, "higher in SLE", "lower in SLE"))
res$label <- factor(res$label, levels = res$label[order(res$odds_ratio)])
cols <- c("higher in SLE" = PAL$sle, "lower in SLE" = PAL$healthy, "not significant" = PAL$muted)
g <- ggplot(res, aes(x = odds_ratio, y = label, colour = direction)) +
  geom_vline(xintercept = 1, colour = PAL$axis, linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.8) +
  geom_point(size = 2.6) +
  scale_x_log10() +
  scale_colour_manual(values = cols, name = NULL) +
  labs(
    title = "Cell-type proportions in lupus vs healthy donors",
    subtitle = sprintf("Odds ratio with 95%% CI, adjusted for %s; %d donors, %d sample x batch units",
                       paste(p2$model$covariates, collapse = " and "), length(unique(meta$donor)), nrow(meta)),
    x = "Odds ratio (log scale; >1 = larger share of cells in SLE)", y = NULL
  ) +
  theme_sle() + theme(legend.position = "top", legend.justification = "left")
fig <- file.path(d$man$figures_dir, "07_proportions_forest.png")
ggsave(fig, g, width = 7, height = 0.32 * nrow(res) + 1.8, dpi = 150, bg = PAL$surface)

# ---- report ---------------------------------------------------------------------------
out <- res[, c("cell_type", "mean_pct_healthy", "mean_pct_SLE", "odds_ratio", "ci_low", "ci_high", "fdr")]
out$fdr <- fmt_p(out$fdr)
names(out) <- c("cell type", "mean % healthy", "mean % SLE", "odds ratio", "95% CI low", "95% CI high", "FDR")
L <- c(
  "# Phase 2: cell-type proportions\n",
  sprintf("- Units: %d sample x batch units from %d donors (%d SLE, %d healthy); cell counts from all cells",
          nrow(meta), length(unique(meta$donor)), length(unique(meta$donor[is_sle])),
          length(unique(meta$donor[!is_sle]))),
  sprintf("- Model: limma on logit(proportion) ~ %s", paste(colnames(design)[-1], collapse = " + ")),
  sprintf("- Within-donor correlation between technical replicates: %s",
          ifelse(is.na(corr), "not needed (no repeated donors)", sprintf("%.2f", corr))),
  "", md_table(out), "",
  "Odds ratio > 1 means that cell type makes up a larger share of a lupus donor's PBMCs.",
  "Proportions are compositional: if one cell type expands, others shrink as a share even when",
  "their absolute numbers do not change, so results are read together, not one by one.",
  "", sprintf("Figure: `%s`", basename(fig))
)
write_md(file.path(d$man$reports_dir, "07_proportions.md"), L)
