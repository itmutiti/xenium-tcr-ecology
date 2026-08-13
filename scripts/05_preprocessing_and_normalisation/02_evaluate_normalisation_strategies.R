#!/usr/bin/env Rscript
# `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` -- 02_evaluate_normalisation_strategies.R
#
# Compares `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`'s four layers (counts, lognorm, pearson_residuals,
# detected) on two axes: replicate stability (pseudobulk correlation
# between the two runs of each of the 7 technical-replicate pairs) and
# negative-control-probe behaviour (Spearman correlation between each
# method's per-cell mean value and that cell's raw control-probe ratio --
# lower magnitude indicates less residual technical noise surviving
# normalisation).
#
# The actual metric computation runs in Python
# (_02_compute_normalization_benchmark_metrics.py, invoked below via
# system2()), not R: this comparison needs direct access to
# analysis_ready.h5ad (1.12M cells x 623 genes x 4 layers), and no R
# HDF5/AnnData reader is available in this project's environment. This script owns the
# actual blueprint-mandated deliverable (the PDF report) and reads back
# the Python helper's output.
#
# This script reports the comparison; it does not itself select the
# method to use -- which layer becomes the primary analysis
# representation for `05_preprocessing_and_normalisation/03_calculate_program_scores.py` onward is recorded as a human decision,
# not decided silently here.
#
# Primary output: reports/preprocess/normalisation_benchmark.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(ggplot2)
  library(tidyr)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

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

METHOD_LABELS <- c(
  counts = "Raw counts (log1p pseudobulk)",
  lognorm = "Log-normalisation",
  pearson_residuals = "Pearson residuals",
  detected = "Binary detection"
)

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  helper_script <- file.path(project_root, "scripts", "05_preprocessing_and_normalisation", "_02_compute_normalization_benchmark_metrics.py")
  cat("[INFO] Running Python metric computation helper...\n")
  exit_code <- system2("python3", c(shQuote(helper_script), "--project-root", shQuote(project_root)))
  if (exit_code != 0) {
    stop("Python metric computation helper failed -- see its output above.")
  }

  replicate_stability <- read_parquet(file.path(project_root, "reports", "preprocess", "normalisation_benchmark_replicate_stability.parquet"))
  technical_noise <- read_parquet(file.path(project_root, "reports", "preprocess", "normalisation_benchmark_technical_noise.parquet"))

  methods <- names(METHOD_LABELS)
  replicate_long <- pivot_longer(
    replicate_stability,
    cols = paste0(methods, "_replicate_r"),
    names_to = "method", values_to = "replicate_r"
  )
  replicate_long$method <- sub("_replicate_r$", "", replicate_long$method)
  replicate_long$method_label <- METHOD_LABELS[replicate_long$method]

  technical_noise$method_label <- METHOD_LABELS[technical_noise$method]
  technical_noise$abs_rho <- abs(technical_noise$spearman_rho_vs_control_probe_ratio)

  output_dir <- file.path(project_root, "reports", "preprocess")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Side-by-side, not stacked: for 2 panels a single moderate-width
  # row keeps the aspect ratio closer to square than either a wide
  # 2-in-a-row or a tall 2x1 stack would.
  open_publication_pdf(file.path(output_dir, "normalisation_benchmark.pdf"), width = 15.5, height = 7.0)

  p1 <- ggplot(replicate_long, aes(x = method_label, y = replicate_r)) +
    geom_boxplot(fill = "white", colour = PUB_COLORS$not_significant, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, colour = PUB_COLORS$primary_analysis, size = 2.2, alpha = 0.8) +
    labs(
      subtitle = "Pearson r between the two runs' pseudobulk\nprofiles (n = 7 replicate pairs)",
      x = NULL, y = "Replicate Pearson r"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  p2 <- ggplot(technical_noise, aes(x = reorder(method_label, abs_rho), y = abs_rho)) +
    geom_col(fill = PUB_COLORS$sensitivity_analysis, width = 0.6) +
    coord_flip() +
    labs(
      subtitle = "|Spearman rho| vs. control-probe ratio\n(lower is better)",
      x = NULL, y = "|Spearman rho|"
    ) +
    theme_publication()

  compose_panels(list(p1, p2), ncol = 2)

  dev.off()

  summary_table <- merge(
    aggregate(replicate_r ~ method, replicate_long, median),
    technical_noise[, c("method", "abs_rho")],
    by = "method"
  )
  colnames(summary_table) <- c("method", "median_replicate_r", "abs_technical_noise_rho")
  write.table(
    summary_table,
    file.path(project_root, "data", "derived", "normalisation_benchmark_summary.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  cat("[OK]   Benchmark summary:\n")
  print(summary_table)
  cat(sprintf(
    "[OK]   Wrote %s and data/derived/normalisation_benchmark_summary.tsv. This script does not select a winner.\n",
    file.path(output_dir, "normalisation_benchmark.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
