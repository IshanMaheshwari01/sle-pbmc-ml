# Shared helpers for the phase-2 R scripts.
# Input comes from scripts/06_export_pseudobulk.py (a folder with counts.mtx.gz,
# profiles.csv, genes.csv and manifest.yaml).

suppressPackageStartupMessages({
  library(Matrix)
  library(yaml)
  library(ggplot2)
})

script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else getwd()
}

# First command-line argument: the phase-2 data folder.
phase2_dir <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1) return(normalizePath(a[1]))
  normalizePath(file.path(script_dir(), "..", "data", "processed", "phase2"))
}

load_phase2 <- function(dir = phase2_dir()) {
  man <- yaml::read_yaml(file.path(dir, "manifest.yaml"))
  if (!file.exists(man$counts)) stop("Counts not found: ", man$counts, ". Run `make export` first.")
  counts <- as(Matrix::readMM(gzfile(man$counts)), "CsparseMatrix")
  profiles <- read.csv(man$profiles, stringsAsFactors = FALSE, check.names = FALSE)
  genes <- read.csv(man$genes, stringsAsFactors = FALSE)
  rownames(counts) <- genes$gene_id
  colnames(counts) <- profiles$profile
  rownames(profiles) <- profiles$profile
  profiles$in_cohort <- as.logical(profiles$in_cohort)
  profiles$low_cell_count <- as.logical(profiles$low_cell_count)
  dir.create(man$reports_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(man$figures_dir, recursive = TRUE, showWarnings = FALSE)
  list(counts = counts, profiles = profiles, genes = genes, man = man, p2 = man$phase2)
}

# Same validated palette as the Python figures.
PAL <- list(
  surface = "#fcfcfb", ink = "#0b0b0b", ink2 = "#52514e", muted = "#898781",
  grid = "#e1e0d9", axis = "#c3c2b7", healthy = "#2a78d6", sle = "#eb6834",
  neutral = "#f0efec", blue_ramp = c("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"),
  red_ramp = c("#f6d5d5", "#eea3a2", "#e34948", "#c03636", "#982a2a")
)

theme_sle <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = PAL$surface, colour = NA),
      panel.background = element_rect(fill = PAL$surface, colour = NA),
      panel.grid.major = element_line(colour = PAL$grid, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(colour = PAL$muted),
      axis.title = element_text(colour = PAL$ink2),
      plot.title = element_text(colour = PAL$ink, size = base_size + 1.5, hjust = 0),
      plot.title.position = "plot",
      plot.subtitle = element_text(colour = PAL$ink2, hjust = 0),
      strip.text = element_text(colour = PAL$ink, hjust = 0),
      legend.title = element_text(colour = PAL$ink2),
      legend.text = element_text(colour = PAL$ink2)
    )
}

slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))

short_ct <- function(x) {
  x <- sub(", alpha-beta T cell", " T cell", x)
  x <- sub("-positive", "+", x)
  x
}

fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-4, "< 1e-4", formatC(p, format = "g", digits = 2)))

md_table <- function(df, digits = 3) {
  df <- as.data.frame(df, check.names = FALSE)
  for (j in seq_along(df)) if (is.numeric(df[[j]])) df[[j]] <- signif(df[[j]], digits)
  head <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(trimws(format(r)), collapse = " | "), " |"))
  c(head, sep, rows)
}

write_md <- function(path, lines) {
  writeLines(lines, path)
  cat(paste(lines, collapse = "\n"), "\n")
  cat("Written:", path, "\n")
}

# Design with disease (healthy as reference) plus covariates that vary in the data.
make_design <- function(meta, covariates) {
  meta$disease_group <- factor(meta$disease_group, levels = c("healthy", "SLE"))
  use <- covariates[vapply(covariates, function(v) length(unique(meta[[v]])) > 1, logical(1))]
  for (v in use) meta[[v]] <- factor(meta[[v]])
  f <- as.formula(paste("~", paste(c("disease_group", use), collapse = " + ")))
  X <- model.matrix(f, data = meta)
  # drop columns that are not estimable (e.g. a batch level with no variation left)
  qx <- qr(X)
  X[, qx$pivot[seq_len(qx$rank)], drop = FALSE]
}
