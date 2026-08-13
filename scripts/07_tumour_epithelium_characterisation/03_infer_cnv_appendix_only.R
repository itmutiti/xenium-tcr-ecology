#!/usr/bin/env Rscript
# `07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R` -- 03_infer_cnv_appendix_only.R
#
# Attempts CNV-like inference from the panel's biological-gene features
# (395 of 399 with resolved genomic coordinates -- see
# scripts/07_tumour_epithelium_characterisation/_03_prepare_cnv_inputs.py); restricted to an appendix
# sensitivity output and barred from any primary claim path, per the
# blueprint's own framing.
#
# Panel-density constraint, checked before choosing a method: the 395
# resolved genes spread across 23 autosomes/X
# (chrY excluded -- only 2 genes, unusable) give a median of ~16
# genes/chromosome (range 4-43). This is far too sparse for sub-
# chromosomal windowed smoothing (the standard inferCNV/CopyKAT approach,
# which assumes hundreds-to-thousands of genes per chromosome) to be
# anything but spurious precision. CNV signal is therefore computed at
# whole-chromosome resolution only -- an explicit, appendix-tier
# limitation, not a hidden one.
#
# Method: for each epithelial cell, relative expression per gene = that
# cell's lognorm value minus its own patient's non-epithelial (immune/
# stromal) reference-cell mean for that gene (patient-matched reference,
# standard CNV-inference practice, avoids a cross-patient technical-batch
# confound). Per-chromosome CNV score = mean relative expression across
# that chromosome's genes. Per-cell CNV burden = variance of the
# per-chromosome scores across chromosomes (an aneuploid cell should show
# more between-chromosome deviation than a diploid one). CNV burden is
# correlated against `07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`'s malignancy_probability as a
# partial, appendix-tier cross-validation -- not treated as an independent
# confirmatory result on its own, consistent with this script's barred-
# from-primary-claims status.
#
# The actual data export runs in Python (_03_prepare_cnv_inputs.py,
# invoked below via system2()): this script cannot read .h5ad directly (no
# R HDF5/AnnData reader available).
#
# Primary output: reports/tumour/cnv_appendix.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

# Real chromosomes only, in conventional karyotype order -- chrY excluded
# per the module-level comment above (2 genes, not usable for a
# per-chromosome mean).
CHROMOSOME_ORDER <- c(as.character(1:22), "X")

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

parse_project_root_arg <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- which(args == "--project-root")
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  NULL
}

# A general (non-Xenium-specific) boolean-as-string interop hazard, already
# documented and fixed once per patient in this project
# -- not needed here since
# no boolean columns cross the Python/R boundary in this script's inputs,
# but `as_bool()` is intentionally not reintroduced speculatively; noted
# only so a future reader does not wonder why it is absent.

genes_by_chromosome <- function(gene_coords) {
  usable <- gene_coords[gene_coords$chromosome %in% CHROMOSOME_ORDER, ]
  split(usable$gene, factor(usable$chromosome, levels = CHROMOSOME_ORDER))
}

compute_relative_expression <- function(epi_expr, reference_baseline, genes) {
  # epi_expr: cells x genes matrix (+ a parallel patient_id vector).
  # reference_baseline: patients x genes matrix, row-indexed by patient_id.
  # Returns a cells x genes matrix of (cell value - that cell's own
  # patient's reference mean), vectorised over cells via row-matched
  # baseline lookup rather than a per-cell loop.
  patient_rows <- match(epi_expr$patient_id, rownames(reference_baseline))
  if (anyNA(patient_rows)) {
    stop("Some cells' patient_id has no matching reference_baseline row.")
  }
  baseline_matrix <- as.matrix(reference_baseline[patient_rows, genes, drop = FALSE])
  cell_matrix <- as.matrix(epi_expr[, genes, drop = FALSE])
  cell_matrix - baseline_matrix
}

compute_chromosome_cnv_scores <- function(relative_expression, chrom_gene_list) {
  # relative_expression: cells x genes matrix. Returns cells x chromosomes
  # matrix of mean relative expression per chromosome.
  scores <- sapply(chrom_gene_list, function(genes) {
    genes_present <- intersect(genes, colnames(relative_expression))
    if (length(genes_present) == 0) return(rep(NA_real_, nrow(relative_expression)))
    if (length(genes_present) == 1) return(relative_expression[, genes_present])
    rowMeans(relative_expression[, genes_present, drop = FALSE])
  })
  matrix(scores, nrow = nrow(relative_expression), dimnames = list(rownames(relative_expression), names(chrom_gene_list)))
}

compute_cnv_burden <- function(chromosome_scores) {
  apply(chromosome_scores, 1, function(row) stats::var(row, na.rm = TRUE))
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  helper_script <- file.path(project_root, "scripts", "07_tumour_epithelium_characterisation", "_03_prepare_cnv_inputs.py")
  cat("[INFO] Running Python CNV-input preparation helper...\n")
  exit_code <- system2("python3", c(shQuote(helper_script), "--project-root", shQuote(project_root)))
  if (exit_code != 0) {
    stop("Python CNV-input preparation helper failed -- see its output above.")
  }

  gene_coords <- read.delim(file.path(project_root, "references", "gene_genomic_coordinates.tsv"), stringsAsFactors = FALSE)
  epi_expr <- as.data.frame(read_parquet(file.path(project_root, "data", "derived", "cnv_epithelial_expression.parquet")))
  reference_baseline <- as.data.frame(read_parquet(file.path(project_root, "data", "derived", "cnv_reference_baseline.parquet")))
  rownames(reference_baseline) <- reference_baseline$patient_id
  reference_baseline$patient_id <- NULL

  chrom_gene_list <- genes_by_chromosome(gene_coords)
  usable_genes <- unlist(chrom_gene_list, use.names = FALSE)
  usable_genes <- intersect(usable_genes, colnames(epi_expr))
  usable_genes <- intersect(usable_genes, colnames(reference_baseline))
  chrom_gene_list <- lapply(chrom_gene_list, function(g) intersect(g, usable_genes))
  chrom_gene_list <- chrom_gene_list[lengths(chrom_gene_list) > 0]

  relative_expression <- compute_relative_expression(epi_expr, reference_baseline, usable_genes)
  chromosome_scores <- compute_chromosome_cnv_scores(relative_expression, chrom_gene_list)
  cnv_burden <- compute_cnv_burden(chromosome_scores)

  cnv_results <- data.frame(
    cell_id = rownames(epi_expr),
    patient_id = epi_expr$patient_id,
    malignancy_probability = epi_expr$malignancy_probability,
    cnv_burden = cnv_burden
  )
  cnv_results <- cbind(cnv_results, as.data.frame(chromosome_scores))

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(cnv_results, file.path(output_dir_data, "cnv_scores.parquet"))

  overall_cor <- suppressWarnings(stats::cor.test(cnv_results$cnv_burden, cnv_results$malignancy_probability, method = "spearman"))

  per_patient_cor <- do.call(rbind, lapply(split(cnv_results, cnv_results$patient_id), function(g) {
    if (nrow(g) < 10 || stats::sd(g$cnv_burden, na.rm = TRUE) == 0) {
      return(data.frame(patient_id = g$patient_id[1], n_cells = nrow(g), spearman_rho = NA_real_, pval = NA_real_))
    }
    test <- suppressWarnings(stats::cor.test(g$cnv_burden, g$malignancy_probability, method = "spearman"))
    data.frame(patient_id = g$patient_id[1], n_cells = nrow(g), spearman_rho = unname(test$estimate), pval = test$p.value)
  }))
  rownames(per_patient_cor) <- NULL

  output_dir_reports <- file.path(project_root, "reports", "tumour")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  open_publication_pdf(file.path(output_dir_reports, "cnv_appendix.pdf"), width = 17, height = 12)

  chrom_long <- data.frame(
    cell_id = rep(cnv_results$cell_id, times = ncol(chromosome_scores)),
    malignancy_probability = rep(cnv_results$malignancy_probability, times = ncol(chromosome_scores)),
    chromosome = rep(colnames(chromosome_scores), each = nrow(chromosome_scores)),
    score = as.vector(chromosome_scores)
  )
  chrom_long$chromosome <- factor(chrom_long$chromosome, levels = colnames(chromosome_scores))
  chrom_long$malignancy_tier <- cut(chrom_long$malignancy_probability, breaks = c(-Inf, 1/3, 2/3, Inf), labels = c("low", "mid", "high"))

  # outlier.shape = NA, not just a small outlier.size: with ~9.1M rows
  # (397,247 cells x 23 chromosomes), ggplot must otherwise individually
  # plot every outlier point, which is the actual bottleneck that made an
  # earlier run of this script exceed a 5-minute timeout -- the boxplot
  # summary statistics themselves are unaffected, only the (uninformative
  # at this density) individual outlier markers are skipped.
  p1 <- ggplot(chrom_long, aes(x = chromosome, y = score, fill = malignancy_tier)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), linewidth = 0.35) +
    scale_fill_manual(values = c("low" = OKABE_ITO$blue, "mid" = OKABE_ITO$grey, "high" = OKABE_ITO$vermillion), name = "Malignancy\ntier") +
    labs(
      subtitle = sprintf(
        "Whole-chromosome resolution only, median %d genes/chromosome (appendix, not a primary claim)",
        stats::median(lengths(chrom_gene_list))
      ),
      x = "Chromosome", y = "Mean relative expression\nvs. patient-matched reference"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  print(p1)

  # geom_bin2d (built into ggplot2, no extra package), not geom_point +
  # geom_smooth(loess): loess scales poorly (effectively O(n^2)) and was
  # the other bottleneck behind the same timeout noted above at
  # n=397,247 points -- a 2D-binned density plot scales to this n trivially
  # and represents an overplotted relationship more accurately at
  # this density than a raw scatter would. The already-computed
  # Spearman rho (robust to nonlinearity) in the subtitle carries the
  # quantitative relationship; no smoother line is added back.
  p2 <- ggplot(cnv_results, aes(x = malignancy_probability, y = cnv_burden)) +
    geom_bin2d(bins = 60) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    labs(
      subtitle = sprintf("Overall Spearman rho = %.3f (p %s)", unname(overall_cor$estimate),
                          ifelse(overall_cor$p.value < 1e-4, "< 1e-4", sprintf("= %.4f", overall_cor$p.value))),
      x = "Malignancy probability", y = "CNV burden\n(variance of per-chromosome scores)"
    ) +
    theme_publication()
  print(p2)

  dev.off()

  cat("[OK]   CNV appendix complete.\n")
  cat(sprintf("[OK]   Overall CNV burden vs malignancy_probability: Spearman rho = %.4f, p %s\n",
              unname(overall_cor$estimate), ifelse(overall_cor$p.value < 1e-4, "< 1e-4", sprintf("= %.4f", overall_cor$p.value))))
  cat("[OK]   Per-patient correlations:\n")
  print(per_patient_cor)
  cat(sprintf("[OK]   Wrote %s, %s\n",
              file.path(output_dir_reports, "cnv_appendix.pdf"),
              file.path(output_dir_data, "cnv_scores.parquet")))
  cat("[OK]   STATUS: APPENDIX ONLY -- whole-chromosome resolution (395 genes, median ~16/chromosome), not sub-chromosomal; not a primary claim path.\n")
}

if (sys.nframe() == 0) {
  main()
}
