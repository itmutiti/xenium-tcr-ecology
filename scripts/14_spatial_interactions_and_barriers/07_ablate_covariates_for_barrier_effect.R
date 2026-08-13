#!/usr/bin/env Rscript
# `14_spatial_interactions_and_barriers/07_ablate_covariates_for_barrier_effect.R`
#
# Decomposes the gap between the prespecified Q3 covariate-adjusted
# suppressive-myeloid barrier effect (`03_model_barrier_topology_by_
# structure.R`: estimate=-0.343, LRT p=0.0069) and the weak, non-
# significant raw bivariate correlation reported by `06_benchmark_
# against_published_barrier_studies.R` (r=-0.079, p=0.332):
# which adjustment covariate(s) account for the difference.
#
# Method: adds each of `03`'s 9 adjustment covariates (4 state-
# composition fractions, 5 niche-archetype fractions) to a barrier-only
# mixed model one at a time, tracking the
# `suppressive_myeloid_barrier_fraction` fixed-effect estimate, SE, and
# Wald p-value at each step; also fits the two covariate blocks
# (state-only, niche-only) and the full model (state+niche, identical
# formula to `03`) to test which covariate class accounts for the
# effect. Reuses `03`'s data-loading, filtering, and niche-composition-
# pivoting logic, so the analysis unit is identical (n=152 clone-sections
# after barrier-reachability exclusion). Added after the pipeline was
# initially complete; see `docs/analysis_amendments.md`.
#
# Primary output: reports/interactions/barrier_covariate_ablation.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(lme4)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

STATE_COVARIATES <- c("cytotoxic_fraction", "exhausted_fraction", "cycling_fraction", "treg_fraction")
BARRIER_COVARIATES <- c("fibroblast_barrier_fraction", "suppressive_myeloid_barrier_fraction")
FOCAL_BARRIER_COVARIATE <- "suppressive_myeloid_barrier_fraction"

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

#' Per-(clone_id, section_id) fractional membership across each archetype
#' in `levels`; identical logic to `03_model_barrier_topology_by_
#' structure.R`'s `compute_niche_composition_wide`, reproduced here by
#' value rather than cross-sourced between sibling scripts.
compute_niche_composition_wide <- function(domain_data, levels) {
  keys <- unique(domain_data[, c("clone_id", "section_id")])
  counts <- table(
    interaction(domain_data$clone_id, domain_data$section_id, drop = TRUE, lex.order = TRUE),
    factor(domain_data$archetype, levels = levels)
  )
  totals <- rowSums(counts)
  fractions <- counts / totals
  key_labels <- interaction(keys$clone_id, keys$section_id, drop = TRUE, lex.order = TRUE)
  fractions_df <- as.data.frame.matrix(fractions)
  colnames(fractions_df) <- paste0("niche_archetype_", levels, "_fraction")
  fractions_df$clone_id <- keys$clone_id[match(rownames(fractions), key_labels)]
  fractions_df$section_id <- keys$section_id[match(rownames(fractions), key_labels)]
  rownames(fractions_df) <- NULL
  fractions_df
}

#' Pure, testable: Nakagawa & Schielzeth (2013) marginal R^2, identical
#' to `03`'s own helper.
compute_marginal_r2 <- function(model) {
  fixed_pred <- predict(model, re.form = NA)
  var_fixed <- stats::var(fixed_pred)
  vc <- as.data.frame(lme4::VarCorr(model))
  var_random_total <- sum(vc$vcov)
  var_fixed / (var_fixed + var_random_total)
}

#' Fits `engagement_ratio ~ BARRIER_COVARIATES + adjustment_covariates +
#' (1 | patient_id/clone_id)` and extracts the focal barrier covariate's
#' fixed-effect estimate, SE, a 95% Wald CI, a normal-approximation Wald
#' p-value, and marginal R^2 -- one row of the ablation table.
#' `lme4::lmer` (unlike `lmerTest`) does not compute degrees-of-freedom-
#' based p-values by default; a standard normal-approximation Wald
#' p-value (`2*pnorm(-|z|)`) is used
#' instead, consistent with the Wald CI already reported and with `03_
#' model_barrier_topology_by_structure.R`'s own choice not to depend on
#' `lmerTest`.
fit_and_extract_barrier_effect <- function(data, adjustment_covariates, step_label) {
  rhs_terms <- c(BARRIER_COVARIATES, adjustment_covariates)
  formula <- stats::as.formula(sprintf("engagement_ratio ~ %s + (1 | patient_id/clone_id)", paste(rhs_terms, collapse = " + ")))
  model <- lme4::lmer(formula, data = data, REML = FALSE)
  coefs <- summary(model)$coefficients
  estimate <- coefs[FOCAL_BARRIER_COVARIATE, "Estimate"]
  se <- coefs[FOCAL_BARRIER_COVARIATE, "Std. Error"]
  data.frame(
    step = step_label,
    n_adjustment_covariates = length(adjustment_covariates),
    estimate = estimate,
    se = se,
    ci_low = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se,
    p_value = 2 * stats::pnorm(-abs(estimate / se)),
    marginal_r2 = compute_marginal_r2(model)
  )
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  engagement_path <- file.path(project_root, "data", "derived", "clone_tumour_engagement.parquet")
  state_path <- file.path(project_root, "data", "derived", "clone_state_composition.parquet")
  barrier_path <- file.path(project_root, "data", "derived", "clone_barrier_metrics.parquet")
  tissue_domains_path <- file.path(project_root, "data", "derived", "tissue_domains.parquet")
  resolved_calls_path <- file.path(project_root, "data", "derived", "tcr_resolved_calls.parquet")
  high_confidence_clones_path <- file.path(project_root, "data", "releases", "v1_tcr_calls", "high_confidence_clones.parquet")
  sample_manifest_path <- file.path(project_root, "metadata", "sample_manifest.tsv")

  for (p in c(engagement_path, state_path, barrier_path, tissue_domains_path, resolved_calls_path, high_confidence_clones_path, sample_manifest_path)) {
    if (!file.exists(p)) stop(sprintf("'%s' not found.", p))
  }

  engagement <- as.data.frame(read_parquet(engagement_path))
  state <- as.data.frame(read_parquet(state_path))
  barrier <- as.data.frame(read_parquet(barrier_path))
  tissue_domains <- as.data.frame(read_parquet(tissue_domains_path))
  resolved <- as.data.frame(read_parquet(resolved_calls_path))
  high_confidence_clones <- as.data.frame(read_parquet(high_confidence_clones_path))
  sample_manifest <- read.delim(sample_manifest_path, sep = "\t", stringsAsFactors = FALSE)

  primary_sections <- sample_manifest$section_id[sample_manifest$included_in_primary_hnscc_cohort == "True"]
  clonal <- resolved[resolved$resolution %in% c("singlet", "low_confidence"), ]
  clonal$clone_id <- clonal$detected_probes
  clonal <- clonal[clonal$clone_id %in% high_confidence_clones$clone_id & clonal$section_id %in% primary_sections, ]

  domain_data <- merge(clonal[, c("cell_id", "clone_id", "section_id")], tissue_domains[, c("cell_id", "section_id", "archetype")], by = c("cell_id", "section_id"))
  archetype_levels <- sort(unique(tissue_domains$archetype))
  niche <- compute_niche_composition_wide(domain_data, archetype_levels)
  niche_covariates <- paste0("niche_archetype_", archetype_levels[-1], "_fraction")

  data <- merge(engagement[, c("clone_id", "section_id", "patient_id", "engagement_ratio")], state[, c("clone_id", "section_id", STATE_COVARIATES)], by = c("clone_id", "section_id"))
  data <- merge(data, barrier[, c("clone_id", "section_id", BARRIER_COVARIATES)], by = c("clone_id", "section_id"))
  data <- merge(data, niche[, c("clone_id", "section_id", niche_covariates)], by = c("clone_id", "section_id"))
  data <- data[stats::complete.cases(data[, c(BARRIER_COVARIATES, niche_covariates)]), ]

  cat(sprintf(
    "[INFO] Analysis unit (identical to `03_model_barrier_topology_by_structure.R`): n=%d clone-section rows, %d clones, %d patients.\n",
    nrow(data), length(unique(data$clone_id)), length(unique(data$patient_id))
  ))

  # Structured ablation: barrier-only (the mixed-model equivalent of the
  # raw correlation, but accounting for patient/clone clustering, unlike
  # the simple Pearson r in `06_benchmark_against_published_barrier_
  # studies.R`), then each of the 9 adjustment covariates added one at a
  # time, then the two covariate blocks alone, then the full model
  # (identical to `03`'s own).
  all_adjustment_covariates <- c(STATE_COVARIATES, niche_covariates)
  steps <- list()
  steps[["barrier_only"]] <- character(0)
  for (cov in all_adjustment_covariates) steps[[paste0("+ ", cov)]] <- cov
  steps[["state_block_only"]] <- STATE_COVARIATES
  steps[["niche_block_only"]] <- niche_covariates
  steps[["full_state_and_niche"]] <- all_adjustment_covariates

  ablation <- do.call(rbind, lapply(names(steps), function(step_label) {
    fit_and_extract_barrier_effect(data, steps[[step_label]], step_label)
  }))
  rownames(ablation) <- NULL

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(ablation, file.path(output_dir_data, "barrier_covariate_ablation.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "interactions")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(output_dir_reports, "barrier_covariate_ablation.pdf"), width = 11, height = 7)

  single_steps <- ablation[!(ablation$step %in% c("barrier_only", "state_block_only", "niche_block_only", "full_state_and_niche")), ]
  single_steps$step <- factor(single_steps$step, levels = single_steps$step)
  p1 <- ggplot(single_steps, aes(x = step, y = estimate)) +
    geom_col(fill = "darkred") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(data = ablation[ablation$step == "barrier_only", ], aes(yintercept = estimate), linetype = "dotted", color = "grey40") +
    labs(
      title = sprintf("Covariate ablation: %s effect, one adjustment covariate added at a time", FOCAL_BARRIER_COVARIATE),
      subtitle = "Dotted line = barrier-only estimate (no adjustment covariates, mixed-model analog of the raw correlation)",
      x = "Single adjustment covariate added to the barrier-only model", y = "Barrier fixed-effect estimate (95% Wald CI)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  print(p1)

  block_steps <- ablation[ablation$step %in% c("barrier_only", "state_block_only", "niche_block_only", "full_state_and_niche"), ]
  block_steps$step <- factor(block_steps$step, levels = c("barrier_only", "state_block_only", "niche_block_only", "full_state_and_niche"))
  p2 <- ggplot(block_steps, aes(x = step, y = estimate)) +
    geom_col(fill = "steelblue") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = sprintf("Block-level ablation: which covariate class unmasks the %s effect?", FOCAL_BARRIER_COVARIATE),
      x = "Adjustment covariate block", y = "Barrier fixed-effect estimate (95% Wald CI)"
    ) +
    theme_minimal()
  print(p2)

  dev.off()

  cat("[INFO] Single-covariate ablation (barrier estimate after adding each covariate alone):\n")
  for (i in seq_len(nrow(single_steps))) {
    cat(sprintf("[INFO]   %-45s %.4f (p=%.4f)\n", single_steps$step[i], single_steps$estimate[i], single_steps$p_value[i]))
  }
  cat("[INFO] Block-level ablation:\n")
  for (i in seq_len(nrow(block_steps))) {
    cat(sprintf("[INFO]   %-25s %.4f (%.4f-%.4f), p=%.4f\n", block_steps$step[i], block_steps$estimate[i], block_steps$ci_low[i], block_steps$ci_high[i], block_steps$p_value[i]))
  }
  cat(sprintf(
    "[OK]   Barrier covariate ablation complete. Wrote %s, %s\n",
    file.path(output_dir_data, "barrier_covariate_ablation.parquet"),
    file.path(output_dir_reports, "barrier_covariate_ablation.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
