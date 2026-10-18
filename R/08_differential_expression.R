#!/usr/bin/env Rscript
# Phase 2, step 3: pseudobulk differential expression, lupus vs healthy, per cell type.
#
# Primary method: limma-voom on pseudobulk counts (TMM-normalised), model
#     expression ~ disease + batch + ancestry
# with duplicateCorrelation blocked on donor, because a sample sequenced in two
# batches gives two correlated profiles. Moderated t-tests, Benjamini-Hochberg FDR.
#
# Robustness check: DESeq2 on one profile per sample (the batch with most cells),
# same covariates. The two methods make different assumptions; genes both call
# are the ones to trust.
#
# Pathways: MSigDB Hallmark sets tested with limma::cameraPR on the moderated
# t-statistics (competitive test that allows for correlation between genes).
#
# Usage: Rscript R/08_differential_expression.R [data/processed/phase2]

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "common.R"))
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})

d <- load_phase2()
p2 <- d$p2
fdr_cut <- p2$de$fdr
lfc_cut <- p2$de$min_abs_logfc
pr <- d$profiles[d$profiles$in_cohort & d$profiles$n_cells >= p2$de$min_cells, ]
sym <- setNames(d$genes$symbol, d$genes$gene_id)
de_dir <- file.path(d$man$reports_dir, "08_de")
dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)

# ---- gene sets ----------------------------------------------------------------------
read_gmt <- function(path) {
  lines <- strsplit(readLines(path), "\t")
  setNames(lapply(lines, function(x) x[-(1:2)]), vapply(lines, `[`, "", 1))
}
load_gene_sets <- function() {
  gmt <- p2$enrichment$gmt
  if (!is.null(gmt) && nzchar(gmt)) {
    if (!grepl("^/", gmt)) gmt <- file.path(script_dir(), "..", gmt)
    if (file.exists(gmt)) return(list(sets = read_gmt(gmt), source = basename(gmt)))
    message("[warn] gene-set file not found: ", gmt)
  }
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    df <- tryCatch(msigdbr::msigdbr(species = "human", collection = p2$enrichment$collection),
                   error = function(e) tryCatch(
                     msigdbr::msigdbr(species = "Homo sapiens", category = p2$enrichment$collection),
                     error = function(e2) NULL))
    if (!is.null(df) && nrow(df)) {
      return(list(sets = split(df$gene_symbol, df$gs_name),
                  source = sprintf("msigdbr %s, collection %s", packageVersion("msigdbr"), p2$enrichment$collection)))
    }
  }
  message("[warn] no gene sets available: install msigdbr or run `make genesets`; skipping pathway tests")
  NULL
}
gs <- load_gene_sets()

# ---- per cell type ----------------------------------------------------------------------
types <- names(sort(table(pr$cell_type), decreasing = TRUE))
summary_rows <- list(); camera_rows <- list(); ifn_rows <- list(); fits <- list()

for (ct in types) {
  m <- pr[pr$cell_type == ct, ]
  n_don <- tapply(m$donor, m$disease_group, function(x) length(unique(x)))
  if (any(is.na(n_don[c("SLE", "healthy")])) || min(n_don[c("SLE", "healthy")]) < p2$de$min_profiles_per_group) {
    message(sprintf("Skipping %s: too few donors (%s)", ct, paste(names(n_don), n_don, collapse = ", ")))
    next
  }
  cat(sprintf("\n== %s: %d profiles, %d SLE / %d healthy donors\n", ct, nrow(m), n_don["SLE"], n_don["healthy"]))
  y <- DGEList(d$counts[, m$profile, drop = FALSE])
  design <- make_design(m, p2$model$covariates)
  keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  coef <- "disease_groupSLE"

  corr <- NA
  if (anyDuplicated(m$donor)) {
    v <- voom(y, design)
    corr <- duplicateCorrelation(v, design, block = m$donor)$consensus.correlation
    v <- voom(y, design, block = m$donor, correlation = corr)
    corr <- duplicateCorrelation(v, design, block = m$donor)$consensus.correlation
    fit <- lmFit(v, design, block = m$donor, correlation = corr)
  } else {
    v <- voom(y, design)
    fit <- lmFit(v, design)
  }
  fit <- eBayes(fit, robust = TRUE)
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
  tt <- data.frame(gene_id = rownames(tt), symbol = sym[rownames(tt)], tt, check.names = FALSE)
  fits[[ct]] <- tt

  # ---- DESeq2 check -------------------------------------------------------------------------
  rho <- NA; overlap <- NA; n_deseq <- NA
  if (isTRUE(p2$de$run_deseq2) && requireNamespace("DESeq2", quietly = TRUE)) {
    one <- m[order(m$sample_id, -m$n_cells), ]
    one <- one[!duplicated(one$sample_id), ]
    res2 <- tryCatch({
      cd <- one[, c("disease_group", p2$model$covariates)]
      cd$disease_group <- factor(cd$disease_group, levels = c("healthy", "SLE"))
      use <- p2$model$covariates[vapply(p2$model$covariates, function(v) length(unique(cd[[v]])) > 1, logical(1))]
      for (vv in use) cd[[vv]] <- factor(cd[[vv]])
      f <- as.formula(paste("~", paste(c(use, "disease_group"), collapse = " + ")))
      cnt <- as.matrix(d$counts[rownames(y), one$profile, drop = FALSE])
      storage.mode(cnt) <- "integer"
      dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(cnt, cd, f))
      dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
      as.data.frame(DESeq2::results(dds, name = "disease_group_SLE_vs_healthy"))
    }, error = function(e) { message("[warn] DESeq2 failed for ", ct, ": ", conditionMessage(e)); NULL })
    if (!is.null(res2)) {
      res2 <- data.frame(gene_id = rownames(res2), symbol = sym[rownames(res2)], res2, check.names = FALSE)
      write.csv(res2, file.path(de_dir, paste0(slug(ct), "_deseq2.csv")), row.names = FALSE)
      shared <- intersect(tt$gene_id, res2$gene_id)
      rho <- suppressWarnings(cor(tt[match(shared, tt$gene_id), "logFC"],
                                  res2[match(shared, res2$gene_id), "log2FoldChange"],
                                  method = "spearman", use = "complete.obs"))
      sig_l <- tt$gene_id[tt$adj.P.Val < fdr_cut]
      sig_d <- res2$gene_id[!is.na(res2$padj) & res2$padj < fdr_cut]
      n_deseq <- length(sig_d)
      overlap <- if (length(sig_l)) 100 * mean(sig_l %in% sig_d) else NA
      tt$deseq2_log2FC <- res2$log2FoldChange[match(tt$gene_id, res2$gene_id)]
      tt$deseq2_padj <- res2$padj[match(tt$gene_id, res2$gene_id)]
      fits[[ct]] <- tt
    }
  }
  write.csv(tt, file.path(de_dir, paste0(slug(ct), "_limma.csv")), row.names = FALSE)

  sig <- tt$adj.P.Val < fdr_cut
  summary_rows[[ct]] <- data.frame(
    cell_type = ct, sle_donors = n_don["SLE"], healthy_donors = n_don["healthy"], profiles = nrow(m),
    genes_tested = nrow(tt), donor_corr = round(corr, 2),
    de_up = sum(sig & tt$logFC > 0), de_down = sum(sig & tt$logFC < 0),
    notable_up = sum(sig & tt$logFC >= lfc_cut), notable_down = sum(sig & tt$logFC <= -lfc_cut),
    deseq2_de = n_deseq, deseq2_rho = round(rho, 2), pct_limma_de_also_deseq2 = round(overlap, 1),
    top_up = paste(head(tt$symbol[sig & tt$logFC > 0], 6), collapse = ", "),
    top_down = paste(head(tt$symbol[sig & tt$logFC < 0], 6), collapse = ", "),
    check.names = FALSE
  )

  ifn <- tt[tt$symbol %in% p2$ifn_genes, c("symbol", "logFC", "adj.P.Val")]
  if (nrow(ifn)) ifn_rows[[ct]] <- data.frame(cell_type = ct, ifn)

  if (!is.null(gs)) {
    idx <- ids2indices(gs$sets, tt$symbol, remove.empty = TRUE)
    idx <- idx[lengths(idx) >= (if (is.null(p2$enrichment$min_set_size)) 10 else p2$enrichment$min_set_size)]
    if (length(idx)) {
      cam <- cameraPR(setNames(tt$t, tt$gene_id), idx)
      cam$FDR <- p.adjust(cam$PValue, "BH")
      camera_rows[[ct]] <- data.frame(cell_type = ct, pathway = rownames(cam), cam, check.names = FALSE)
    }
  }
}

if (!length(summary_rows)) stop("No cell type had enough donors in both groups.")
summ <- do.call(rbind, summary_rows)
write.csv(summ, file.path(d$man$reports_dir, "08_de_summary.csv"), row.names = FALSE)
cam_all <- if (length(camera_rows)) do.call(rbind, camera_rows) else NULL
if (!is.null(cam_all)) write.csv(cam_all, file.path(d$man$reports_dir, "08_pathways_camera.csv"), row.names = FALSE)

# ---- figures ------------------------------------------------------------------------------
figs <- character()
# 1. number of DE genes per cell type, up and down
bars <- rbind(
  data.frame(cell_type = summ$cell_type, n = summ$notable_up, dir = "higher in SLE"),
  data.frame(cell_type = summ$cell_type, n = -summ$notable_down, dir = "lower in SLE")
)
bars$label <- factor(short_ct(bars$cell_type), levels = rev(short_ct(summ$cell_type)))
g1 <- ggplot(bars, aes(x = n, y = label, fill = dir)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0, colour = PAL$axis) +
  geom_text(data = bars[bars$n != 0, ], aes(label = abs(n), hjust = ifelse(n > 0, -0.25, 1.25)),
            size = 2.8, colour = PAL$ink2) +
  scale_fill_manual(values = c("higher in SLE" = PAL$sle, "lower in SLE" = PAL$healthy), name = NULL) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.15)) +
  labs(title = "Differentially expressed genes per cell type (lupus vs healthy)",
       subtitle = sprintf("limma-voom, FDR < %.2f and |log2 fold change| >= %.1f", fdr_cut, lfc_cut),
       x = "Genes", y = NULL) +
  theme_sle() + theme(legend.position = "top", legend.justification = "left")
f1 <- file.path(d$man$figures_dir, "08_de_counts.png")
ggsave(f1, g1, width = 7, height = 0.33 * nrow(summ) + 1.8, dpi = 150, bg = PAL$surface)
figs <- c(figs, f1)

# 2. volcano for the cell type with the most DE genes
top_ct <- summ$cell_type[which.max(summ$de_up + summ$de_down)]
tv <- fits[[top_ct]]
tv$cls <- ifelse(tv$adj.P.Val < fdr_cut & tv$logFC >= lfc_cut, "higher in SLE",
                 ifelse(tv$adj.P.Val < fdr_cut & tv$logFC <= -lfc_cut, "lower in SLE", "other"))
lab <- head(tv[tv$cls != "other", ], 12)
g2 <- ggplot(tv, aes(logFC, -log10(P.Value), colour = cls)) +
  geom_point(size = 0.9, alpha = 0.7) +
  geom_text(data = lab, aes(label = symbol), size = 2.6, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = c("higher in SLE" = PAL$sle, "lower in SLE" = PAL$healthy, "other" = PAL$axis), name = NULL) +
  labs(title = sprintf("%s: lupus vs healthy", short_ct(top_ct)),
       subtitle = "Pseudobulk limma-voom, adjusted for batch and ancestry; labels = strongest hits",
       x = "log2 fold change (SLE / healthy)", y = "-log10 p-value") +
  theme_sle() + theme(legend.position = "top", legend.justification = "left")
f2 <- file.path(d$man$figures_dir, paste0("08_volcano_", slug(top_ct), ".png"))
ggsave(f2, g2, width = 6.5, height = 5, dpi = 150, bg = PAL$surface)
figs <- c(figs, f2)

# 3. limma vs DESeq2 agreement for the same cell type
if (!is.null(tv$deseq2_log2FC)) {
  g3 <- ggplot(tv[!is.na(tv$deseq2_log2FC), ], aes(logFC, deseq2_log2FC)) +
    geom_abline(slope = 1, intercept = 0, colour = PAL$axis) +
    geom_point(size = 0.7, alpha = 0.5, colour = PAL$healthy) +
    labs(title = sprintf("%s: two methods, same answer?", short_ct(top_ct)),
         subtitle = sprintf("Spearman rho = %.2f across %d genes",
                            summ$deseq2_rho[summ$cell_type == top_ct], sum(!is.na(tv$deseq2_log2FC))),
         x = "limma-voom log2 fold change", y = "DESeq2 log2 fold change") +
    coord_equal() + theme_sle()
  f3 <- file.path(d$man$figures_dir, "08_limma_vs_deseq2.png")
  ggsave(f3, g3, width = 5, height = 5, dpi = 150, bg = PAL$surface)
  figs <- c(figs, f3)
}

# 4. pathways x cell types heatmap (signed -log10 FDR)
if (!is.null(cam_all)) {
  cam_all$score <- ifelse(cam_all$Direction == "Up", 1, -1) * -log10(pmax(cam_all$FDR, 1e-20))
  hit <- unique(cam_all$pathway[cam_all$FDR < p2$enrichment$fdr])
  if (length(hit)) {
    best <- tapply(abs(cam_all$score), cam_all$pathway, max)
    hit <- head(names(sort(best[hit], decreasing = TRUE)), 20)
    hm <- cam_all[cam_all$pathway %in% hit, ]
    hm$pathway_label <- factor(sub("^HALLMARK_", "", hm$pathway), levels = rev(sub("^HALLMARK_", "", hit)))
    hm$ct_label <- factor(short_ct(hm$cell_type), levels = short_ct(summ$cell_type))
    lim <- max(abs(hm$score))
    g4 <- ggplot(hm, aes(ct_label, pathway_label, fill = score)) +
      geom_tile(colour = PAL$surface, linewidth = 0.6) +
      geom_text(aes(label = ifelse(FDR < p2$enrichment$fdr, "*", "")), colour = PAL$ink, size = 3.5, vjust = 0.75) +
      scale_fill_gradient2(low = PAL$healthy, mid = PAL$neutral, high = "#e34948", midpoint = 0,
                           limits = c(-lim, lim), name = "signed\n-log10 FDR") +
      labs(title = "Hallmark pathways shifted in lupus, by cell type",
           subtitle = sprintf("cameraPR on limma t-statistics; red = up in SLE, blue = down; * FDR < %.2f",
                              p2$enrichment$fdr), x = NULL, y = NULL) +
      theme_sle() + theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
    f4 <- file.path(d$man$figures_dir, "08_pathways_heatmap.png")
    ggsave(f4, g4, width = 0.55 * length(unique(hm$cell_type)) + 4.5, height = 0.28 * length(hit) + 2.2,
           dpi = 150, bg = PAL$surface)
    figs <- c(figs, f4)
  }
}

# ---- report --------------------------------------------------------------------------------
st <- summ[, c("cell_type", "sle_donors", "healthy_donors", "genes_tested", "de_up", "de_down",
               "notable_up", "notable_down", "deseq2_rho", "pct_limma_de_also_deseq2")]
names(st) <- c("cell type", "SLE donors", "healthy donors", "genes tested", "DE up", "DE down",
               sprintf("up, |lfc|>=%.1f", lfc_cut), sprintf("down, |lfc|>=%.1f", lfc_cut),
               "limma vs DESeq2 rho", "% limma DE also DESeq2 DE")
L <- c(
  "# Phase 2: pseudobulk differential expression\n",
  sprintf("- Cohort: %d donors, profiles with >= %d cells; one sample per donor, technical replicates across batches kept",
          length(unique(pr$donor)), p2$de$min_cells),
  sprintf("- limma-voom ~ disease + %s, duplicateCorrelation on donor; robust eBayes; BH FDR < %.2f",
          paste(p2$model$covariates, collapse = " + "), fdr_cut),
  sprintf("- DESeq2 check: %s", ifelse(isTRUE(p2$de$run_deseq2), "one profile per sample, same covariates", "off")),
  "", md_table(st), "",
  "## Strongest genes per cell type\n",
  md_table(summ[, c("cell_type", "top_up", "top_down")]), ""
)
if (length(ifn_rows)) {
  ifn <- do.call(rbind, ifn_rows)
  wide <- reshape(ifn[, c("cell_type", "symbol", "logFC")], idvar = "cell_type", timevar = "symbol", direction = "wide")
  names(wide) <- sub("^logFC\\.", "", names(wide))
  L <- c(L, "## Interferon genes (log2 fold change, SLE vs healthy)\n",
         "The six genes of the phase-1 interferon score, now tested properly at the donor level.\n",
         md_table(wide, digits = 2), "")
}
if (!is.null(cam_all)) {
  topc <- cam_all[cam_all$FDR < p2$enrichment$fdr, ]
  topc <- topc[order(topc$FDR), ][seq_len(min(25, nrow(topc))), c("cell_type", "pathway", "NGenes", "Direction", "FDR")]
  topc$FDR <- fmt_p(topc$FDR)
  L <- c(L, sprintf("## Pathways (%s), FDR < %.2f\n", gs$source, p2$enrichment$fdr),
         if (nrow(topc)) md_table(topc) else "_None significant._", "")
}
L <- c(L, "## Files\n", sprintf("- Per-gene tables: `08_de/<cell type>_limma.csv` (and `_deseq2.csv`)"),
       sprintf("- Figures: %s", paste0("`", basename(figs), "`", collapse = ", ")))
write_md(file.path(d$man$reports_dir, "08_differential_expression.md"), L)
