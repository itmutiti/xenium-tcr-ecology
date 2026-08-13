#!/usr/bin/env Rscript
# `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R` -- 01_fit_variance_partition_models.R
#
# The prespecified `q2_variance_partition_confirmatory` analysis
# (governance/analysis_registry.tsv): "Clonal spatial/phenotypic
# descriptor variance decomposes into identity, context, and patient
# components with non-trivial, estimable magnitude for each -- not that
# spatial context is the only or dominant source."
#
# Target variable: `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s frozen `ecological_structure_
# score` (data/releases/v1_clone_structure/, taxonomy_version=
# v1_provisional) -- the same continuous score this project's whole
# clone-ecology arc has been building toward, not a freshly chosen
# descriptor. Per-(clone, section) observations (n=261, 158 distinct
# clones, 10 patients).
#
# Three-level variance decomposition (nested random-intercept model,
# `lme4`, the same package already validated for exactly this purpose
# in `10_niche_and_ecosystem_discovery/06_test_patient_recurrence.R`):
#   ecological_structure_score ~ 1 + (1 | patient_id/clone_id)
# -- patient-level random-intercept variance = patient component;
# clone-level random-intercept variance (nested within patient, valid
# because `clone_id` values are never shared across patients by
# construction, `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`, `08_tcr_clonal_analysis/07_build_clone_metadata_table.py`'s patient-restriction invariant) = identity
# component; residual variance (section-to-section variation for the
# same clone) = context component -- a clone's score changing across
# its different spatial locations is exactly what "local
# microenvironmental context" means.
#
# Primary output: reports/clone_ecology/variance_partition_models.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(lme4)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

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

#' Fits the three-level nested variance-partition model.
fit_variance_partition_model <- function(data) {
  lme4::lmer(ecological_structure_score ~ 1 + (1 | patient_id/clone_id), data = data, REML = TRUE)
}

#' Pure, testable: extracts the three named variance components from a
#' fitted `fit_variance_partition_model` object.
extract_variance_partition <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  list(
    patient_var = vc$vcov[vc$grp == "patient_id"],
    identity_var = vc$vcov[vc$grp == "clone_id:patient_id"],
    context_var = vc$vcov[vc$grp == "Residual"]
  )
}

#' Pure, testable: each component's share of total variance.
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

  structure_path <- file.path(project_root, "data", "releases", "v1_clone_structure", "clone_ecological_structure.parquet")
  if (!file.exists(structure_path)) stop(sprintf("'%s' not found. Run `13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py` first.", structure_path))

  data <- as.data.frame(read_parquet(structure_path))
  cat(sprintf("[INFO] Data: n=%d clone-section rows, %d clones, %d patients.\n", nrow(data), length(unique(data$clone_id)), length(unique(data$patient_id))))

  model <- fit_variance_partition_model(data)
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
  # bootMer returns a row of NA for any replicate whose refit failed to
  # converge or errored; na.rm = TRUE below silently excludes those from
  # the CI, so the number of replicates actually informing the CI is
  # tracked and reported explicitly rather than assumed to equal nsim.
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
  write_parquet(result, file.path(output_dir_data, "variance_partition_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "clone_ecology")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(output_dir_reports, "variance_partition_models.pdf"), width = 9, height = 6)

  result$component <- factor(result$component, levels = c("patient", "identity", "context"))
  p1 <- ggplot(result, aes(x = component, y = proportion)) +
    geom_col(fill = "steelblue") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
    labs(
      title = "Q2: clone ecological-structure-score variance partition",
      subtitle = sprintf("n=%d clone-sections, %d clones, %d patients; 95%% bootstrap CIs (%d/%d reps successful)", nrow(data), length(unique(data$clone_id)), length(unique(data$patient_id)), n_boot_success, N_BOOTSTRAP),
      x = "Variance component", y = "Proportion of total variance"
    ) +
    ylim(0, 1) +
    theme_minimal()
  print(p1)

  dev.off()

  cat("[INFO] Variance partition (proportion of total, 95% bootstrap CI):\n")
  for (i in seq_len(nrow(result))) {
    cat(sprintf("[INFO]   %-10s %.3f (%.3f-%.3f)\n", result$component[i], result$proportion[i], result$ci_low[i], result$ci_high[i]))
  }
  cat(sprintf(
    "[OK]   Variance partition complete (%d/%d bootstrap reps successful). Wrote %s, %s\n",
    n_boot_success, N_BOOTSTRAP,
    file.path(output_dir_data, "variance_partition_results.parquet"),
    file.path(output_dir_reports, "variance_partition_models.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
