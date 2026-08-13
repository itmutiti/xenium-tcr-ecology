#!/usr/bin/env Rscript
# `16_external_validation_and_generalisation/03_validate_cell_state_signatures.R` -- 03_validate_cell_state_signatures.R
#
# Tests the `cell_state_signature_generalisation` claim (governance/
# validation_plan.tsv, `16_external_validation_and_generalisation/00_define_validation_claims.py`): whether this project's T-cell
# programme gene signatures (`05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s `PROGRAM_GENE_SETS`) are
# coherently co-expressed modules in GSE139324 (Cillo et al. 2020,
# Immunity) -- a second independent HNSCC scRNA-seq reference,
# deliberately different from GSE103322/Puram et al. 2017 already used
# in `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` (see data/external/scrna/GSE139324/README.md).
#
# Scope limitation: GSE139324 is a CD45+ immune-sorted cohort (no
# stromal/epithelial cells present, by its experimental design) -- this
# milestone can only test the T-cell programmes (cytotoxicity,
# exhaustion, proliferation, treg -- the same 4 programmes
# `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` already tested on GSE103322, reused here for direct
# comparability, per governance/validation_plan.tsv's stated success
# criterion), not the "stromal states" this milestone's own scaffold
# title also names -- no independent stromal reference was acquired in
# `16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`
# (only bulk RNA-seq, not single-cell, which cannot resolve stromal
# cell-level module coherence the same way).
#
# Method: identifies T cells within each of GSE139324's 26 TIL samples
# via canonical CD3D/CD3E/CD3G marker positivity, reads 10x MTX triplets
# directly in R (base `Matrix`/`data.table`, no Python intermediate),
# computes a CPM-like log-normalised expression value per marker gene,
# and tests module coherence (mean pairwise correlation among a
# programme's genes vs. a random background gene pool) -- the same
# methodological question `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`'s `program_transfer.py` asked
# on GSE103322, reimplemented natively in R here (this project's
# convention: R scripts do not source Python modules).
#
# Primary output: reports/validation/signature_validation.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"
N_PERMUTATIONS <- 999
N_BACKGROUND_GENES <- 300
MIN_EXPRESSED_FRACTION <- 0.05
T_CELL_MARKERS <- c("CD3D", "CD3E", "CD3G")

# Exact match to `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s `PROGRAM_GENE_SETS` + Cell Type Annotation's
# `TREG_MARKERS` -- the same 4 programmes `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` tested on
# GSE103322 (`PROGRAMS_TESTED` + Treg tested separately), redefined
# here per this project's R-script-independence convention.
PROGRAM_GENE_SETS <- list(
  cytotoxicity = c("GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "KLRD1", "KLRB1", "KLRC1", "FGFBP2"),
  exhaustion = c("PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "ENTPD1", "TOX"),
  proliferation = c("MKI67", "TOP2A", "PCNA", "CCNB2", "CDK1", "UBE2C", "CENPF"),
  treg = c("FOXP3", "IL2RA", "CTLA4")
)

find_project_root <- function(cli_arg = NULL) {
  check_marker <- function(candidate) file.exists(file.path(candidate, PROJECT_ROOT_MARKER))

  if (!is.null(cli_arg) && nzchar(cli_arg)) {
    candidate <- normalizePath(cli_arg, mustWork = FALSE)
    if (check_marker(candidate)) return(candidate)
    stop(sprintf("Project root supplied via --project-root ('%s') does not contain '%s'.", candidate, PROJECT_ROOT_MARKER))
  }
  env_root <- Sys.getenv(ENV_ROOT_VAR, unset = NA)
  if (!is.na(env_root) && nzchar(env_root)) {
    candidate <- normalizePath(env_root, mustWork = FALSE)
    if (check_marker(candidate)) return(candidate)
    stop(sprintf("Project root supplied via $%s ('%s') does not contain '%s'.", ENV_ROOT_VAR, candidate, PROJECT_ROOT_MARKER))
  }
  here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
  parent <- here
  repeat {
    if (check_marker(parent)) return(parent)
    next_parent <- dirname(parent)
    if (next_parent == parent) break
    parent <- next_parent
  }
  stop(sprintf("Could not locate the project root. None of --project-root, $%s, or a '%s' marker file found by walking up from this script identified it.", ENV_ROOT_VAR, PROJECT_ROOT_MARKER))
}

RNG_SEED <- local({
  # centrally configured (config/config.yaml's default_seed), resolved via
  # a self-contained getwd()-based walk (not find_project_root(), which
  # depends on Rscript's own --file= commandArgs and is unavailable when
  # this script is merely source()-d, e.g. by tests/r_unit/*.R)
  root <- getwd()
  repeat {
    if (file.exists(file.path(root, PROJECT_ROOT_MARKER))) break
    parent <- dirname(root)
    if (parent == root) stop("Could not locate the project root from getwd() to read config/config.yaml for RNG_SEED.")
    root <- parent
  }
  yaml::read_yaml(file.path(root, "config", "config.yaml"))$default_seed
})

parse_project_root_arg <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- which(args == "--project-root")
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  NULL
}

#' Per-cell T-cell flag -- TRUE if any T-cell marker gene (rows of
#' `marker_counts`) has a positive count in that cell (column).
identify_t_cells_from_marker_counts <- function(marker_counts) {
  as.vector(Matrix::colSums(marker_counts > 0) > 0)
}

#' CPM-like normalisation (counts / per-cell total library size x
#' median total, avoiding raw-count correlations being dominated by
#' cell-depth artefacts) then `log1p`.
normalize_counts_cpm_log <- function(counts, total_counts, median_total) {
  log1p(counts / total_counts * median_total)
}

#' Mean of the upper-triangle Pearson correlation matrix among the
#' columns (genes) of a cells x genes normalised expression matrix.
compute_mean_pairwise_correlation <- function(expr_matrix) {
  corr <- stats::cor(expr_matrix, method = "pearson")
  mean(corr[upper.tri(corr)])
}

#' Permutation-based empirical p-value -- draws `n_permutations` random
#' gene sets of size `n_genes` from the columns of `background_expr`,
#' computes each one's mean pairwise correlation, and compares to
#' `observed_correlation`.
compute_module_coherence_pvalue <- function(observed_correlation, background_expr, n_genes, n_permutations = N_PERMUTATIONS) {
  n_background_genes <- ncol(background_expr)
  null_stats <- vapply(seq_len(n_permutations), function(i) {
    idx <- sample.int(n_background_genes, n_genes)
    compute_mean_pairwise_correlation(background_expr[, idx, drop = FALSE])
  }, numeric(1))
  pvalue <- (sum(null_stats >= observed_correlation) + 1) / (n_permutations + 1)
  list(pvalue = pvalue, null_mean = mean(null_stats))
}

#' Reads one 10x CellRanger v2 MTX triplet (genes x cells).
read_10x_triplet <- function(sample_dir, prefix) {
  matrix_path <- file.path(sample_dir, paste0(prefix, "_matrix.mtx.gz"))
  genes_path <- file.path(sample_dir, paste0(prefix, "_genes.tsv.gz"))
  barcodes_path <- file.path(sample_dir, paste0(prefix, "_barcodes.tsv.gz"))
  mat <- Matrix::readMM(gzfile(matrix_path))
  # data.table::fread's direct-gz-read path requires the 'R.utils'
  # package, not installed in this environment; decompressing via a
  # gzfile() connection first avoids the dependency.
  genes <- data.table::fread(text = readLines(gzfile(genes_path)), header = FALSE, col.names = c("ensembl_id", "symbol"))
  barcodes <- data.table::fread(text = readLines(gzfile(barcodes_path)), header = FALSE, col.names = c("barcode"))
  rownames(mat) <- genes$symbol
  colnames(mat) <- barcodes$barcode
  mat
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))
  raw_dir <- file.path(project_root, "data", "external", "scrna", "GSE139324", "raw")
  if (!dir.exists(raw_dir)) stop(sprintf("'%s' not found. Run `16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py` first.", raw_dir))

  all_genes_needed <- unique(unlist(PROGRAM_GENE_SETS))
  til_files <- list.files(raw_dir, pattern = "HNSCC_[0-9]+_TIL_matrix\\.mtx\\.gz$")
  til_prefixes <- sub("_matrix\\.mtx\\.gz$", "", til_files)
  cat(sprintf("[INFO] %d TIL samples found.\n", length(til_prefixes)))

  set.seed(RNG_SEED)
  needed_expr_list <- list()
  background_expr_list <- list()
  background_gene_pool <- NULL

  for (prefix in til_prefixes) {
    mat <- read_10x_triplet(raw_dir, prefix)
    total_counts <- Matrix::colSums(mat)
    marker_counts <- mat[T_CELL_MARKERS[T_CELL_MARKERS %in% rownames(mat)], , drop = FALSE]
    is_t_cell <- identify_t_cells_from_marker_counts(marker_counts)
    if (sum(is_t_cell) == 0) next

    t_cell_mat <- mat[, is_t_cell, drop = FALSE]
    t_cell_totals <- total_counts[is_t_cell]
    median_total <- stats::median(t_cell_totals)

    needed_counts <- as.matrix(t_cell_mat[all_genes_needed, , drop = FALSE])
    needed_norm <- t(apply(needed_counts, 1, normalize_counts_cpm_log, total_counts = t_cell_totals, median_total = median_total))
    needed_expr_list[[prefix]] <- t(needed_norm)

    if (is.null(background_gene_pool)) {
      # Sampled once (first sample), reused across all samples for a
      # consistent background gene set -- a fraction-expressed filter
      # (matching `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`'s fix for the same issue: a random
      # background gene pool can include never-expressed genes, which
      # would bias the null toward zero correlation) applied on this
      # first sample as a proxy.
      candidate_genes <- setdiff(rownames(mat), all_genes_needed)
      expressed_fraction <- Matrix::rowMeans(mat[candidate_genes, is_t_cell, drop = FALSE] > 0)
      eligible <- candidate_genes[expressed_fraction >= MIN_EXPRESSED_FRACTION]
      background_gene_pool <- sample(eligible, min(N_BACKGROUND_GENES, length(eligible)))
    }
    background_counts <- as.matrix(t_cell_mat[background_gene_pool, , drop = FALSE])
    background_norm <- t(apply(background_counts, 1, normalize_counts_cpm_log, total_counts = t_cell_totals, median_total = median_total))
    background_expr_list[[prefix]] <- t(background_norm)
  }

  needed_expr <- do.call(rbind, needed_expr_list)
  background_expr <- do.call(rbind, background_expr_list)
  cat(sprintf("[INFO] Combined T-cell pool: n=%d cells across %d samples.\n", nrow(needed_expr), length(needed_expr_list)))

  rows <- list()
  for (program in names(PROGRAM_GENE_SETS)) {
    genes <- PROGRAM_GENE_SETS[[program]]
    observed_corr <- compute_mean_pairwise_correlation(needed_expr[, genes, drop = FALSE])
    coherence <- compute_module_coherence_pvalue(observed_corr, background_expr, n_genes = length(genes))
    rows[[program]] <- data.frame(
      program = program,
      n_genes = length(genes),
      observed_correlation = observed_corr,
      background_mean_correlation = coherence$null_mean,
      pvalue = coherence$pvalue
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  output_dir_data <- file.path(project_root, "data", "derived")
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(result, file.path(output_dir_data, "gse139324_signature_validation_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "validation")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Narrower canvas so in-image text shrinks less when auto-scaled to
  # the page's ~6.5in text width.
  open_publication_pdf(file.path(output_dir_reports, "signature_validation.pdf"), width = 8.6, height = 7.4)

  plot_data <- result
  plot_data$program <- factor(plot_data$program, levels = plot_data$program[order(plot_data$observed_correlation, decreasing = TRUE)])
  p1 <- ggplot(plot_data, aes(x = program)) +
    geom_col(aes(y = observed_correlation, fill = "Observed"), width = 0.6) +
    geom_point(aes(y = background_mean_correlation, colour = "Background (permuted)"), size = 3.8) +
    scale_fill_manual(values = c("Observed" = PUB_COLORS$primary_analysis), name = NULL) +
    scale_colour_manual(values = c("Background (permuted)" = OKABE_ITO$black), name = NULL) +
    labs(
      subtitle = sprintf("Module coherence in GSE139324 (Cillo et al. 2020)\nn = %d permutations", N_PERMUTATIONS),
      x = NULL, y = "Mean pairwise Pearson correlation"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")
  print(p1)

  dev.off()

  cat("[INFO] Module coherence (GSE139324, second independent reference):\n")
  for (i in seq_len(nrow(result))) {
    cat(sprintf("[INFO]   %-16s r=%.3f (background=%.3f) p=%.4f\n", result$program[i], result$observed_correlation[i], result$background_mean_correlation[i], result$pvalue[i]))
  }
  n_significant <- sum(result$pvalue < 0.05)
  cat(sprintf(
    "[OK]   Signature validation complete: %d/%d programme(s) significant (p<0.05). Wrote %s, %s\n",
    n_significant, nrow(result),
    file.path(output_dir_data, "gse139324_signature_validation_results.parquet"),
    file.path(output_dir_reports, "signature_validation.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
