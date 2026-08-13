#!/usr/bin/env Rscript
# `13_clone_ecology_confirmatory_models/03_test_clone_size_as_confounder.R` -- 03_test_clone_size_as_confounder.R
#
# Tests whether `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`'s variance-partition conclusion survives
# adjustment for clone size (`n_cells`) -- the direct proxy this script
# uses for both "clone size" and "detection power" (a larger clone is,
# by construction, both bigger and more reliably/precisely detected;
# TCR-probe-level detection efficiency itself was already characterised
# separately at the probe level in `08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`'s empirical
# false-positive-rate estimation, a different sensitivity
# dimension from this clone-descriptor-level check, out of this
# milestone's more direct scope).
#
# Method: refits `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`'s model with `log(n_cells + 1)` added
# as a fixed-effect covariate, and compares the random-effect
# variance-component proportions before and after adjustment -- if they
# change only modestly, the original three-component conclusion is robust
# to clone-size confounding; a substantial shift would mean size was
# doing more work than the unadjusted model suggested.
#
# Primary output: reports/clone_ecology/clone_size_sensitivity.pdf

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

#' Size-adjusted variance-partition model: same nested random structure
#' as `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`, plus a fixed effect for log(n_cells+1).
fit_size_adjusted_model <- function(data) {
  data$log_n_cells <- log(data$n_cells + 1)
  lme4::lmer(ecological_structure_score ~ 1 + log_n_cells + (1 | patient_id/clone_id), data = data, REML = TRUE)
}

#' Pure, testable: same variance-component extraction as Phase 13.01
#' (reused by name, redefined here to keep this script self-contained
#' per this project's established no-cross-sourcing convention).
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

#' Pure, testable: absolute shift in each component's proportion
#' between the unadjusted and size-adjusted models.
compute_proportion_shift <- function(unadjusted, adjusted) {
  list(
    patient_shift = adjusted$patient_proportion - unadjusted$patient_proportion,
    identity_shift = adjusted$identity_proportion - unadjusted$identity_proportion,
    context_shift = adjusted$context_proportion - unadjusted$context_proportion
  )
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  structure_path <- file.path(project_root, "data", "releases", "v1_clone_structure", "clone_ecological_structure.parquet")
  if (!file.exists(structure_path)) stop(sprintf("'%s' not found. Run `13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py` first.", structure_path))

  data <- as.data.frame(read_parquet(structure_path))

  unadjusted_model <- lme4::lmer(ecological_structure_score ~ 1 + (1 | patient_id/clone_id), data = data, REML = TRUE)
  unadjusted_vc <- extract_variance_partition(unadjusted_model)
  unadjusted_prop <- compute_variance_proportions(unadjusted_vc$patient_var, unadjusted_vc$identity_var, unadjusted_vc$context_var)

  adjusted_model <- fit_size_adjusted_model(data)
  adjusted_vc <- extract_variance_partition(adjusted_model)
  adjusted_prop <- compute_variance_proportions(adjusted_vc$patient_var, adjusted_vc$identity_var, adjusted_vc$context_var)

  shift <- compute_proportion_shift(unadjusted_prop, adjusted_prop)

  size_fixed_effect <- lme4::fixef(adjusted_model)["log_n_cells"]
  size_fixed_effect_se <- sqrt(diag(vcov(adjusted_model)))["log_n_cells"]

  set.seed(RNG_SEED)
  correlation_test <- stats::cor.test(data$n_cells, data$ecological_structure_score, method = "spearman")

  result <- data.frame(
    component = c("patient", "identity", "context"),
    unadjusted_proportion = c(unadjusted_prop$patient_proportion, unadjusted_prop$identity_proportion, unadjusted_prop$context_proportion),
    size_adjusted_proportion = c(adjusted_prop$patient_proportion, adjusted_prop$identity_proportion, adjusted_prop$context_proportion),
    shift = c(shift$patient_shift, shift$identity_shift, shift$context_shift)
  )

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "clone_size_sensitivity_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "clone_ecology")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Side-by-side, not stacked: for 2 panels a single moderate-width
  # row keeps the aspect ratio closer to square than either a wide
  # 2-in-a-row or a tall 2x1 stack would.
  open_publication_pdf(file.path(output_dir_reports, "clone_size_sensitivity.pdf"), width = 15.5, height = 7.0)

  plot_data <- rbind(
    data.frame(component = result$component, proportion = result$unadjusted_proportion, model = "Unadjusted"),
    data.frame(component = result$component, proportion = result$size_adjusted_proportion, model = "Size-adjusted")
  )
  plot_data$component <- factor(plot_data$component, levels = c("patient", "identity", "context"), labels = c("Patient", "Clonal identity", "Spatial context"))
  p1 <- ggplot(plot_data, aes(x = component, y = proportion, fill = model)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("Unadjusted" = PUB_COLORS$not_significant, "Size-adjusted" = PUB_COLORS$primary_analysis)) +
    labs(
      subtitle = sprintf("log(n_cells+1) coefficient = %.4f (SE %.4f)", size_fixed_effect, size_fixed_effect_se),
      x = NULL, y = "Proportion of total variance", fill = NULL
    ) +
    ylim(0, 1) +
    theme_publication()

  p2 <- ggplot(data, aes(x = n_cells, y = ecological_structure_score)) +
    geom_point(alpha = 0.6, colour = PUB_COLORS$not_significant, size = 2.2) +
    scale_x_log10() +
    geom_smooth(method = "lm", se = TRUE, colour = PUB_COLORS$sensitivity_analysis, linewidth = 1.3) +
    labs(
      subtitle = sprintf("Spearman rho = %.3f, p = %.3f", correlation_test$estimate, correlation_test$p.value),
      x = "Clone size, n cells (log scale)", y = "Ecological-structure score"
    ) +
    theme_publication()

  compose_panels(list(p1, p2), ncol = 2)

  dev.off()

  cat("[INFO] Variance proportion shift (size-adjusted minus unadjusted):\n")
  for (i in seq_len(nrow(result))) {
    cat(sprintf("[INFO]   %-10s %.3f -> %.3f (shift %+.3f)\n", result$component[i], result$unadjusted_proportion[i], result$size_adjusted_proportion[i], result$shift[i]))
  }
  cat(sprintf("[INFO] log(n_cells+1) fixed effect: %.4f (SE %.4f)\n", size_fixed_effect, size_fixed_effect_se))
  cat(sprintf("[INFO] Spearman(n_cells, ecological_structure_score): rho=%.3f, p=%.3f\n", correlation_test$estimate, correlation_test$p.value))
  cat(sprintf(
    "[OK]   Clone-size sensitivity check complete. Wrote %s, %s\n",
    file.path(output_dir_data, "clone_size_sensitivity_results.parquet"),
    file.path(output_dir_reports, "clone_size_sensitivity.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
