#!/usr/bin/env Rscript
# `13_clone_ecology_confirmatory_models/07_test_structure_sensitivity_excluding_cycling.R`
#
# Tests how much Q2's headline variance-partition finding (patient 29.3%,
# identity 20.2%, context 50.4%) depends on `cycling_fraction`, the
# dominant loading (-0.998) on the frozen continuous
# `ecological_structure_score` (`11_clone_spatial_descriptors/06_
# discover_provisional_structure.R`) and the input state with the
# weakest external corroboration of those tested
# (40.3% of this cohort's T cells vs. 14.4% in an independent
# reference).
#
# Scope: a sensitivity check, not a replacement. Does not modify the
# frozen `taxonomy_version=v1_provisional` release, the primary
# `ecological_structure_score`, or any published result. Recomputes the
# same 1-factor model (`06_discover_provisional_structure.R`'s
# `compute_continuous_factor_scores`/`extract_factor_loadings`) on a
# feature set with `cycling_fraction` removed (10 of the original 11
# `PRIMARY_FEATURES`), refits the same variance-partition model
# (`13_clone_ecology_confirmatory_models/01_fit_variance_partition_
# models.R`'s `fit_variance_partition_model`), and reports the
# comparison. Added after the pipeline was initially complete; see
# `docs/analysis_amendments.md`.
#
# Primary output: data/derived/structure_sensitivity_excluding_cycling.parquet,
#                  data/derived/variance_partition_sensitivity_excluding_cycling.parquet

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(lme4)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

# Identical to `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s
# `PRIMARY_FEATURES`, minus `cycling_fraction` -- the feature this
# sensitivity check tests the removal of.
PRIMARY_FEATURES_EXCLUDING_CYCLING <- c(
  "cytotoxic_fraction", "exhausted_fraction", "treg_fraction", "cd4_fraction", "cd8_fraction",
  "shannon_entropy_bits", "engagement_ratio", "dc_engagement_ratio", "macrophage_engagement_ratio",
  "antigen_presentation_score_excess"
)
N_BOOTSTRAP <- 500

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

#' Identical to `06_discover_provisional_structure.R`'s own
#' `compute_continuous_factor_scores` -- reused by value, not imported
#' (see that script's own header for why R scripts in this project do
#' not cross-source one another).
compute_continuous_factor_scores <- function(data, n_factors = 1) {
  data_scaled <- scale(data)
  fa <- factanal(data_scaled, factors = n_factors, scores = "regression")
  as.numeric(fa$scores[, 1])
}

extract_factor_loadings <- function(data, n_factors = 1) {
  data_scaled <- scale(data)
  fa <- factanal(data_scaled, factors = n_factors, scores = "regression")
  loadings_vec <- as.numeric(fa$loadings[, 1])
  data.frame(feature = colnames(data), loading = loadings_vec)[order(-abs(loadings_vec)), ]
}

fit_variance_partition_model <- function(data) {
  lme4::lmer(ecological_structure_score ~ 1 + (1 | patient_id/clone_id), data = data, REML = TRUE)
}

extract_variance_partition <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  list(
    patient_var = vc$vcov[vc$grp == "patient_id"],
    identity_var = vc$vcov[vc$grp == "clone_id:patient_id"],
    context_var = vc$vcov[vc$grp == "Residual"]
  )
}

compute_variance_proportions <- function(patient_var, identity_var, context_var) {
  total <- patient_var + identity_var + context_var
  list(
    patient_proportion = patient_var / total,
    identity_proportion = identity_var / total,
    context_proportion = context_var / total
  )
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  paths <- list(
    spatial = file.path(project_root, "data", "derived", "clone_spatial_descriptors.parquet"),
    state = file.path(project_root, "data", "derived", "clone_state_composition.parquet"),
    engagement = file.path(project_root, "data", "derived", "clone_tumour_engagement.parquet"),
    apc = file.path(project_root, "data", "derived", "clone_apc_support.parquet"),
    barrier = file.path(project_root, "data", "derived", "clone_barrier_metrics.parquet"),
    original_variance_partition = file.path(project_root, "data", "derived", "variance_partition_results.parquet")
  )
  for (p in paths) {
    if (!file.exists(p)) stop(sprintf("'%s' not found. Run `11_clone_spatial_descriptors/06_discover_provisional_structure.R` and `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R` first.", p))
  }

  spatial <- as.data.frame(read_parquet(paths$spatial))
  state <- as.data.frame(read_parquet(paths$state))
  engagement <- as.data.frame(read_parquet(paths$engagement))
  apc <- as.data.frame(read_parquet(paths$apc))
  barrier <- as.data.frame(read_parquet(paths$barrier))
  original <- as.data.frame(read_parquet(paths$original_variance_partition))

  key <- c("clone_id", "section_id", "patient_id", "n_cells")
  merged <- Reduce(function(x, y) merge(x, y, by = key), list(spatial, state, engagement, apc, barrier))

  primary_data <- merged[stats::complete.cases(merged[, PRIMARY_FEATURES_EXCLUDING_CYCLING]), ]
  primary_matrix <- as.matrix(primary_data[, PRIMARY_FEATURES_EXCLUDING_CYCLING])

  set.seed(RNG_SEED)
  primary_data$ecological_structure_score <- compute_continuous_factor_scores(primary_matrix, n_factors = 1)
  loadings <- extract_factor_loadings(primary_matrix, n_factors = 1)

  cat(sprintf("[INFO] Analysis unit (cycling_fraction excluded): n=%d clone-section rows, %d clones, %d patients.\n", nrow(primary_data), length(unique(primary_data$clone_id)), length(unique(primary_data$patient_id))))
  cat("[INFO] Top factor loadings without cycling_fraction (|loading| descending):\n")
  for (i in seq_len(min(5, nrow(loadings)))) {
    cat(sprintf("[INFO]   %-35s %+.3f\n", loadings$feature[i], loadings$loading[i]))
  }

  model <- fit_variance_partition_model(primary_data)
  vc <- extract_variance_partition(model)
  proportions <- compute_variance_proportions(vc$patient_var, vc$identity_var, vc$context_var)

  set.seed(RNG_SEED)
  boot <- lme4::bootMer(
    model,
    FUN = function(m) {
      vc_b <- extract_variance_partition(m)
      p_b <- compute_variance_proportions(vc_b$patient_var, vc_b$identity_var, vc_b$context_var)
      c(patient = p_b$patient_proportion, identity = p_b$identity_proportion, context = p_b$context_proportion)
    },
    nsim = N_BOOTSTRAP, use.u = FALSE, type = "parametric"
  )
  # Same convergence-tracking as 01_fit_variance_partition_models.R's identical
  # bootMer call -- found missing here during the second Vast.ai clean-room run:
  # `original`'s schema (read from that script's output) includes these two
  # columns, so rbind()-ing it against this script's `result` for the
  # comparison plot below failed with a column-count mismatch whenever this
  # script's own result lacked them.
  n_boot_success <- sum(stats::complete.cases(boot$t))
  if (n_boot_success < N_BOOTSTRAP) {
    cat(sprintf("[WARN] %d/%d bootstrap replicate(s) failed to converge and were excluded from the CI.\n", N_BOOTSTRAP - n_boot_success, N_BOOTSTRAP))
  }
  ci <- apply(boot$t, 2, function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE))

  result <- data.frame(
    component = c("patient", "identity", "context"),
    variance = c(vc$patient_var, vc$identity_var, vc$context_var),
    proportion = c(proportions$patient_proportion, proportions$identity_proportion, proportions$context_proportion),
    ci_low = ci[1, ],
    ci_high = ci[2, ],
    n_bootstrap_requested = N_BOOTSTRAP,
    n_bootstrap_successful = n_boot_success
  )

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(primary_data[, c("clone_id", "section_id", "patient_id", "n_cells", "ecological_structure_score")], file.path(output_dir_data, "structure_sensitivity_excluding_cycling.parquet"))
  write_parquet(loadings, file.path(output_dir_data, "structure_sensitivity_excluding_cycling_loadings.parquet"))
  write_parquet(result, file.path(output_dir_data, "variance_partition_sensitivity_excluding_cycling.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "clone_ecology")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(output_dir_reports, "structure_sensitivity_excluding_cycling.pdf"), width = 9, height = 6)

  loadings$feature <- factor(loadings$feature, levels = loadings$feature)
  p1 <- ggplot(loadings, aes(x = loading, y = feature, fill = loading > 0)) +
    geom_col() +
    geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
    labs(
      title = "Sensitivity check: continuous structure axis with cycling_fraction excluded",
      subtitle = sprintf("1-factor model, n=%d clone-sections, %d features (cycling_fraction removed)", nrow(primary_data), length(PRIMARY_FEATURES_EXCLUDING_CYCLING)),
      x = "Loading on the continuous structure axis", y = "Feature"
    ) +
    theme_minimal()
  print(p1)

  comparison <- rbind(
    cbind(original, analysis = "original (11 features, incl. cycling_fraction)"),
    cbind(result, analysis = "sensitivity (10 features, cycling_fraction excluded)")
  )
  comparison$component <- factor(comparison$component, levels = c("patient", "identity", "context"))
  p2 <- ggplot(comparison, aes(x = component, y = proportion, fill = analysis)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), width = 0.2) +
    labs(
      title = "Q2 variance partition: original vs. cycling_fraction-excluded sensitivity check",
      subtitle = "95% bootstrap CIs (500 reps), both analyses",
      x = "Variance component", y = "Proportion of total variance", fill = NULL
    ) +
    ylim(0, 1) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p2)

  dev.off()

  cat("[INFO] Variance partition without cycling_fraction (proportion of total, 95% bootstrap CI):\n")
  for (i in seq_len(nrow(result))) {
    cat(sprintf("[INFO]   %-10s %.3f (%.3f-%.3f)  [original: %.3f (%.3f-%.3f)]\n", result$component[i], result$proportion[i], result$ci_low[i], result$ci_high[i], original$proportion[i], original$ci_low[i], original$ci_high[i]))
  }
  cat(sprintf(
    "[OK]   Sensitivity check complete. Wrote %s, %s, %s, %s\n",
    file.path(output_dir_data, "structure_sensitivity_excluding_cycling.parquet"),
    file.path(output_dir_data, "structure_sensitivity_excluding_cycling_loadings.parquet"),
    file.path(output_dir_data, "variance_partition_sensitivity_excluding_cycling.parquet"),
    file.path(output_dir_reports, "structure_sensitivity_excluding_cycling.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
