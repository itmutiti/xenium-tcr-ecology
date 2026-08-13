#!/usr/bin/env Rscript
# `11_clone_spatial_descriptors/06_discover_provisional_structure.R` -- 06_discover_provisional_structure.R
#
# Fits whichever structure `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R` supports, without imposing
# categories the data don't show. `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s result: both the primary
# (n=261) and sensitivity (n=120) analyses concluded continuous -- the
# prespecified null (a continuous latent axis) was retained, discrete
# taxonomy was not supported by either the BIC or gap-statistic
# criterion together. This script therefore fits a continuous 1-factor
# latent score, not a discrete cluster taxonomy -- fitting discrete
# clusters here would directly contradict `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s prespecified
# conclusion.
#
# Reuses `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s feature set (`PRIMARY_FEATURES`,
# `SENSITIVITY_EXTRA_FEATURES`, identical values, kept in sync by hand)
# -- the continuous score reported here must be built from the same
# feature set already validated as the better fit in `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`, not a
# freshly re-derived one. Not sourced directly from 05's script file:
# `find_project_root()`'s path-detection relies on `commandArgs()`,
# which reflects the calling process's arguments, not the sourced file's
# -- cross-sourcing would silently break when this script is itself
# `source()`d from a testthat test file rather than run via `Rscript`
# directly. Every other R script in this project independently redefines
# this same small boilerplate for exactly this reason (consistent with
# `08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`, 9.x, 10.02, 10.06: no existing R script in this
# project cross-sources another script's file).
#
# Primary output: data/derived/clone_ecological_structure.parquet

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

PRIMARY_FEATURES <- c(
  "cytotoxic_fraction", "exhausted_fraction", "cycling_fraction", "treg_fraction", "cd4_fraction", "cd8_fraction",
  "shannon_entropy_bits", "engagement_ratio", "dc_engagement_ratio", "macrophage_engagement_ratio",
  "antigen_presentation_score_excess"
)
SENSITIVITY_EXTRA_FEATURES <- c(
  "dispersion_rarefied", "clark_evans_index_rarefied", "domain_richness_rarefied",
  "fibroblast_barrier_fraction", "vascular_barrier_fraction", "suppressive_myeloid_barrier_fraction"
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

#' 1-factor continuous latent score per observation, via
#' `factanal(..., scores="regression")` (Thomson's regression method,
#' base R). Pure, testable: returns a numeric vector of factor scores,
#' same length/order as `nrow(data)`.
compute_continuous_factor_scores <- function(data, n_factors = 1) {
  data_scaled <- scale(data)
  fa <- factanal(data_scaled, factors = n_factors, scores = "regression")
  as.numeric(fa$scores[, 1])
}

#' Loadings table for a fitted 1-factor model -- which features drive
#' the continuous axis, and in which direction. Pure, testable.
extract_factor_loadings <- function(data, n_factors = 1) {
  data_scaled <- scale(data)
  fa <- factanal(data_scaled, factors = n_factors, scores = "regression")
  loadings_vec <- as.numeric(fa$loadings[, 1])
  data.frame(feature = colnames(data), loading = loadings_vec)[order(-abs(loadings_vec)), ]
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  paths <- list(
    spatial = file.path(project_root, "data", "derived", "clone_spatial_descriptors.parquet"),
    state = file.path(project_root, "data", "derived", "clone_state_composition.parquet"),
    engagement = file.path(project_root, "data", "derived", "clone_tumour_engagement.parquet"),
    apc = file.path(project_root, "data", "derived", "clone_apc_support.parquet"),
    barrier = file.path(project_root, "data", "derived", "clone_barrier_metrics.parquet"),
    structure_test = file.path(project_root, "data", "derived", "clone_structure_test_results.parquet")
  )
  for (p in paths) {
    if (!file.exists(p)) stop(sprintf("'%s' not found. Run `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`-11.05 first.", p))
  }

  structure_test <- as.data.frame(read_parquet(paths$structure_test))
  if (!all(structure_test$final_decision == "continuous")) {
    stop(
      "`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s result is not uniformly 'continuous' -- this script only implements the continuous branch. ",
      "If `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R` has been re-run and now supports discrete structure, this script must be extended (or a discrete ",
      "counterpart written) before proceeding -- fitting continuous structure here would silently contradict a ",
      "discrete-structure finding, which is not a decision this script may make unilaterally."
    )
  }

  spatial <- as.data.frame(read_parquet(paths$spatial))
  state <- as.data.frame(read_parquet(paths$state))
  engagement <- as.data.frame(read_parquet(paths$engagement))
  apc <- as.data.frame(read_parquet(paths$apc))
  barrier <- as.data.frame(read_parquet(paths$barrier))

  key <- c("clone_id", "section_id", "patient_id", "n_cells")
  merged <- Reduce(function(x, y) merge(x, y, by = key), list(spatial, state, engagement, apc, barrier))

  primary_data <- merged[stats::complete.cases(merged[, PRIMARY_FEATURES]), ]
  primary_matrix <- as.matrix(primary_data[, PRIMARY_FEATURES])

  set.seed(RNG_SEED)
  primary_data$ecological_structure_score <- compute_continuous_factor_scores(primary_matrix, n_factors = 1)
  loadings <- extract_factor_loadings(primary_matrix, n_factors = 1)

  # Robustness check (not a fresh, independent analysis): does the
  # primary factor score correlate with a score refit on the richer
  # sensitivity feature set (`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s second, complete-case
  # population), for the clone-sections present in both?
  sensitivity_features <- c(PRIMARY_FEATURES, SENSITIVITY_EXTRA_FEATURES)
  sensitivity_data <- merged[stats::complete.cases(merged[, sensitivity_features]), ]
  set.seed(RNG_SEED)
  sensitivity_data$ecological_structure_score_sensitivity <- compute_continuous_factor_scores(as.matrix(sensitivity_data[, sensitivity_features]), n_factors = 1)

  comparison <- merge(
    primary_data[, c("clone_id", "section_id", "ecological_structure_score")],
    sensitivity_data[, c("clone_id", "section_id", "ecological_structure_score_sensitivity")],
    by = c("clone_id", "section_id")
  )
  robustness_correlation <- if (nrow(comparison) >= 3) stats::cor(comparison$ecological_structure_score, comparison$ecological_structure_score_sensitivity, method = "spearman") else NA_real_

  result <- primary_data[, c("clone_id", "section_id", "patient_id", "n_cells", "ecological_structure_score")]

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "clone_ecological_structure.parquet"))
  write_parquet(loadings, file.path(output_dir_data, "clone_ecological_structure_loadings.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "clone_ecology")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Side-by-side, not stacked: for 2 panels a single moderate-width
  # row keeps the aspect ratio closer to square than either a wide
  # 2-in-a-row or a tall 2x1 stack would. A gets more width than B
  # (11 long feature-name labels vs. a simple histogram).
  open_publication_pdf(file.path(output_dir_reports, "provisional_structure.pdf"), width = 15.5, height = 7.5)

  loadings$feature <- factor(loadings$feature, levels = loadings$feature)
  p1 <- ggplot(loadings, aes(x = loading, y = feature, fill = loading > 0)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    scale_fill_manual(values = c("FALSE" = OKABE_ITO$blue, "TRUE" = OKABE_ITO$vermillion), guide = "none") +
    labs(
      subtitle = sprintf("1-factor model, n = %d clone-sections, %d features", nrow(primary_data), length(PRIMARY_FEATURES)),
      x = "Loading on the continuous structure axis", y = NULL
    ) +
    theme_publication()

  p2 <- ggplot(primary_data, aes(x = ecological_structure_score)) +
    geom_histogram(bins = 30, fill = PUB_COLORS$primary_analysis, colour = "white", linewidth = 0.3) +
    labs(
      subtitle = sprintf(
        "n = %d clone-sections\nSpearman vs. sensitivity-feature refit: %s",
        nrow(primary_data), ifelse(is.na(robustness_correlation), "NA", sprintf("%.3f (n = %d)", robustness_correlation, nrow(comparison)))
      ),
      x = "Ecological-structure score", y = "Clone-sections"
    ) +
    theme_publication()

  compose_panels(list(p1, p2), ncol = 2, widths = c(1.3, 1))

  dev.off()

  cat("[INFO] Top factor loadings (|loading| descending):\n")
  for (i in seq_len(min(5, nrow(loadings)))) {
    cat(sprintf("[INFO]   %-35s %+.3f\n", loadings$feature[i], loadings$loading[i]))
  }
  cat(sprintf(
    "[OK]   %d clone-section(s) scored on the continuous ecological structure axis. Primary/sensitivity score correlation: %s. Wrote %s, %s, %s\n",
    nrow(result),
    ifelse(is.na(robustness_correlation), "NA", sprintf("%.3f (n=%d)", robustness_correlation, nrow(comparison))),
    file.path(output_dir_data, "clone_ecological_structure.parquet"),
    file.path(output_dir_data, "clone_ecological_structure_loadings.parquet"),
    file.path(output_dir_reports, "provisional_structure.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
