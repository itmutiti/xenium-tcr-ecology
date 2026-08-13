#!/usr/bin/env Rscript
# `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R` -- 03_model_barrier_topology_by_structure.R
#
# The prespecified `q3_barrier_topology_confirmatory` analysis
# (governance/analysis_registry.tsv): "Fibroblast/myeloid barrier
# topology ... explains residual variance in clone-tumour engagement
# beyond cell-intrinsic state and niche composition alone." Unit of
# analysis: clone, nested in patient and section -- the same
# (clone_id, section_id) unit already used throughout Clone Spatial Descriptors/13
# (n=261 clone-section rows before any barrier-reachability exclusion,
# not an artefact).
#
# Registered prerequisite ("requires a calibrated degree-preserving
# null (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`) to have passed for the barrier-graph null model
# specifically"): `11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`'s `clone_barrier_metrics.parquet` already
# reuses `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s calibrated null-model type (constrained
# label-permutation within a section, N_PERMUTATIONS=199) -- see
# src/xenium_tcr_ecology/clone_ecology/barrier_metrics.py's module
# docstring. This script consumes that already-computed, already-
# calibrated per-clone barrier data directly; it does not recompute or
# re-calibrate the null.
#
# Target: `engagement_ratio` (data/derived/clone_tumour_
# engagement.parquet, `11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`). Covariates: "cell-intrinsic
# state" = `11_clone_spatial_descriptors/01_compute_clone_cell_state_composition.py`'s per-clone T-cell state-composition
# fractions (data/derived/clone_state_composition.parquet); "niche
# composition" = per-clone fractional membership across Phase 10.02/
# 10.03's discovered k=6 neighbourhood-archetype tissue domains
# (data/derived/tissue_domains.parquet), computed here for the first
# time at clone level (no prior phase needed a per-clone niche
# composition, only per-clone niche diversity --
# `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`'s `domain_richness_rarefied`). Predictor of
# interest: "barrier topology" = `11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`'s
# `fibroblast_barrier_fraction` and `suppressive_myeloid_barrier_
# fraction` (vascular deliberately excluded -- the registered
# hypothesis names fibroblast/myeloid specifically, not vascular).
#
# Method: nested fixed-effects likelihood-ratio test (REML=FALSE,
# standard for LRT on fixed effects), baseline (state + niche) vs full
# (state + niche + barrier), same nested random-intercept structure
# already validated in `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R` (`1 | patient_id/clone_id`). This
# 2-df joint LRT (both barrier covariates together) answers the
# registered hypothesis as written (analysis_registry.tsv); a second,
# separate 1-df LRT isolates suppressive_myeloid_barrier_fraction alone
# (state + niche + fibroblast vs. + myeloid) since that narrower claim,
# not the joint block, is what release gatekeeping and the paper's
# headline result test (`17_statistical_closure_and_release/01_control_multiplicity_and_report_effects.R`
# Bonferroni-corrects this 1-df LRT's p-value, not the joint one). Both
# LRTs are persisted (`barrier_topology_model_results.parquet`'s
# `lrt_chisq`/`lrt_df`/`lrt_pvalue`, myeloid-only, and the joint
# `comparison` values are reported in this script's own log/figure but
# not currently persisted, unchanged from before). Effect size: Nakagawa
# & Schielzeth (2013) marginal R^2 (fixed-effect variance / total
# variance) for baseline and full models, and their difference -- this
# project has no `MuMIn` package available, so marginal R^2 is computed
# directly from `VarCorr` and `predict(..., re.form=NA)` rather than via
# a dependency.
#
# Primary output: reports/interactions/barrier_topology_models.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(lme4)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

N_BOOTSTRAP <- 500
STATE_COVARIATES <- c("cytotoxic_fraction", "exhausted_fraction", "cycling_fraction", "treg_fraction")
BARRIER_COVARIATES <- c("fibroblast_barrier_fraction", "suppressive_myeloid_barrier_fraction")

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

#' Pure, testable: per-(clone_id, section_id) fractional membership
#' across every archetype in `levels`, wide-pivoted as
#' `niche_archetype_<k>_fraction` columns summing to 1.0 per row.
#' `domain_data` has one row per (clone, section)-restricted cell:
#' columns `clone_id`, `section_id`, `archetype`.
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

#' Pure, testable: Nakagawa & Schielzeth (2013) marginal R^2 (fixed-
#' effect variance / total variance) for a fitted `lmer` model with only
#' random-intercept terms (no random slopes) -- this script's model
#' structure.
compute_marginal_r2 <- function(model) {
  fixed_pred <- predict(model, re.form = NA)
  var_fixed <- stats::var(fixed_pred)
  vc <- as.data.frame(lme4::VarCorr(model))
  var_random_total <- sum(vc$vcov)
  var_fixed / (var_fixed + var_random_total)
}

#' Pure, testable: nested-model likelihood-ratio test between a
#' `baseline` and `full` model (both REML=FALSE, fixed-effect nesting
#' required).
compare_nested_models <- function(baseline, full) {
  comparison <- stats::anova(baseline, full)
  list(chisq = comparison$Chisq[2], df = comparison$Df[2], pvalue = comparison$`Pr(>Chisq)`[2])
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

  # Clonal cell membership uses the same restriction `11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`'s
  # barrier_metrics.py already applies (singlet/low_confidence
  # resolution, release clone_ids, primary-cohort sections only) --
  # reproduced here so the niche-composition join uses the identical
  # cell universe as every other `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R` covariate table.
  # `read.delim` reads this Python-written TSV boolean column as the
  # literal strings "True"/"False", not an R logical -- indexing with it
  # silently produces an all-NA section list (a character index vector
  # is treated as a name lookup, not a mask) rather than erroring, so
  # this must be compared explicitly, not relied on to auto-coerce.
  primary_sections <- sample_manifest$section_id[sample_manifest$included_in_primary_hnscc_cohort == "True"]
  clonal <- resolved[resolved$resolution %in% c("singlet", "low_confidence"), ]
  clonal$clone_id <- clonal$detected_probes
  clonal <- clonal[clonal$clone_id %in% high_confidence_clones$clone_id & clonal$section_id %in% primary_sections, ]

  domain_data <- merge(clonal[, c("cell_id", "clone_id", "section_id")], tissue_domains[, c("cell_id", "section_id", "archetype")], by = c("cell_id", "section_id"))
  archetype_levels <- sort(unique(tissue_domains$archetype))
  niche <- compute_niche_composition_wide(domain_data, archetype_levels)
  # Drop the first archetype level as the reference category --
  # fractions sum to 1.0 per row, so including all K would be perfectly
  # collinear with the model intercept (the standard dummy-variable
  # trap).
  niche_covariates <- paste0("niche_archetype_", archetype_levels[-1], "_fraction")

  data <- merge(engagement[, c("clone_id", "section_id", "patient_id", "engagement_ratio")], state[, c("clone_id", "section_id", STATE_COVARIATES)], by = c("clone_id", "section_id"))
  data <- merge(data, barrier[, c("clone_id", "section_id", BARRIER_COVARIATES)], by = c("clone_id", "section_id"))
  data <- merge(data, niche[, c("clone_id", "section_id", niche_covariates)], by = c("clone_id", "section_id"))

  n_before_barrier_exclusion <- nrow(data)
  data <- data[stats::complete.cases(data[, c(BARRIER_COVARIATES, niche_covariates)]), ]
  n_excluded_undefined_barrier <- n_before_barrier_exclusion - nrow(data)

  cat(sprintf(
    "[INFO] Data: n=%d clone-section rows (%d excluded, undefined barrier fraction -- all cells directly tumour-adjacent, zero-length path), %d clones, %d patients.\n",
    nrow(data), n_excluded_undefined_barrier, length(unique(data$clone_id)), length(unique(data$patient_id))
  ))

  baseline_formula <- stats::as.formula(sprintf(
    "engagement_ratio ~ %s + (1 | patient_id/clone_id)",
    paste(c(STATE_COVARIATES, niche_covariates), collapse = " + ")
  ))
  full_formula <- stats::as.formula(sprintf(
    "engagement_ratio ~ %s + (1 | patient_id/clone_id)",
    paste(c(STATE_COVARIATES, niche_covariates, BARRIER_COVARIATES), collapse = " + ")
  ))

  baseline_model <- lme4::lmer(baseline_formula, data = data, REML = FALSE)
  full_model <- lme4::lmer(full_formula, data = data, REML = FALSE)

  # Joint 2-df test of the registered hypothesis (analysis_registry.tsv's
  # q3_barrier_topology_confirmatory: "Fibroblast/myeloid barrier
  # topology... explains residual variance", both covariates together).
  comparison <- compare_nested_models(baseline_model, full_model)
  r2_baseline <- compute_marginal_r2(baseline_model)
  r2_full <- compute_marginal_r2(full_model)

  # Separate 1-df LRT isolating suppressive_myeloid_barrier_fraction
  # specifically (state + niche + fibroblast vs. + myeloid): the joint
  # 2-df test above answers "does the barrier block help," not "is the
  # myeloid coefficient itself significant" -- the latter is what the
  # paper's headline claim (and the Wald test previously used for
  # release gatekeeping) is actually about, so it needs its own,
  # properly nested LRT rather than reusing the joint one.
  myeloid_baseline_formula <- stats::as.formula(sprintf(
    "engagement_ratio ~ %s + (1 | patient_id/clone_id)",
    paste(c(STATE_COVARIATES, niche_covariates, "fibroblast_barrier_fraction"), collapse = " + ")
  ))
  myeloid_baseline_model <- lme4::lmer(myeloid_baseline_formula, data = data, REML = FALSE)
  myeloid_lrt <- compare_nested_models(myeloid_baseline_model, full_model)

  set.seed(RNG_SEED)
  boot <- lme4::bootMer(
    full_model,
    FUN = function(m) lme4::fixef(m)[BARRIER_COVARIATES],
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
  colnames(ci) <- BARRIER_COVARIATES

  full_summary <- summary(full_model)$coefficients
  # lrt_* columns are only populated for suppressive_myeloid_barrier_fraction
  # (NA for fibroblast_barrier_fraction) -- no equivalent single-covariate
  # LRT was computed for fibroblast, since it is not the claim gated at
  # release; see myeloid_lrt above.
  is_myeloid <- BARRIER_COVARIATES == "suppressive_myeloid_barrier_fraction"
  barrier_effects <- data.frame(
    covariate = BARRIER_COVARIATES,
    estimate = full_summary[BARRIER_COVARIATES, "Estimate"],
    se = full_summary[BARRIER_COVARIATES, "Std. Error"],
    ci_low = ci[1, ],
    ci_high = ci[2, ],
    n_bootstrap_requested = N_BOOTSTRAP,
    n_bootstrap_successful = n_boot_success,
    lrt_chisq = ifelse(is_myeloid, myeloid_lrt$chisq, NA_real_),
    lrt_df = ifelse(is_myeloid, myeloid_lrt$df, NA_integer_),
    lrt_pvalue = ifelse(is_myeloid, myeloid_lrt$pvalue, NA_real_)
  )

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(barrier_effects, file.path(output_dir_data, "barrier_topology_model_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "interactions")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(output_dir_reports, "barrier_topology_models.pdf"), width = 9, height = 6)

  plot_data <- barrier_effects
  plot_data$covariate <- factor(plot_data$covariate, levels = BARRIER_COVARIATES)
  p1 <- ggplot(plot_data, aes(x = covariate, y = estimate)) +
    geom_col(fill = "darkred") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = "Q3: barrier-topology fixed effects on clone-tumour engagement",
      subtitle = sprintf(
        "n=%d clone-sections, %d clones, %d patients; barrier-block LRT chisq=%.2f (df=%d), p=%.4f; myeloid-only LRT chisq=%.2f (df=%d), p=%.4f; marginal R^2 baseline=%.3f, full=%.3f; bootstrap CI from %d/%d reps",
        nrow(data), length(unique(data$clone_id)), length(unique(data$patient_id)),
        comparison$chisq, comparison$df, comparison$pvalue,
        myeloid_lrt$chisq, myeloid_lrt$df, myeloid_lrt$pvalue,
        r2_baseline, r2_full, n_boot_success, N_BOOTSTRAP
      ),
      x = "Barrier covariate", y = "Fixed-effect estimate (95% bootstrap CI)"
    ) +
    theme_minimal()
  print(p1)

  dev.off()

  cat(sprintf("[INFO] Nested LRT, barrier block (fibroblast+myeloid beyond state+niche): chisq=%.3f, df=%d, p=%.4f\n", comparison$chisq, comparison$df, comparison$pvalue))
  cat(sprintf("[INFO] Nested LRT, suppressive_myeloid_barrier_fraction only (beyond state+niche+fibroblast): chisq=%.3f, df=%d, p=%.4f\n", myeloid_lrt$chisq, myeloid_lrt$df, myeloid_lrt$pvalue))
  cat(sprintf("[INFO] Marginal R^2: baseline=%.4f, full=%.4f, delta=%.4f\n", r2_baseline, r2_full, r2_full - r2_baseline))
  for (i in seq_len(nrow(barrier_effects))) {
    cat(sprintf(
      "[INFO]   %-38s %.4f (%.4f-%.4f)\n",
      barrier_effects$covariate[i], barrier_effects$estimate[i], barrier_effects$ci_low[i], barrier_effects$ci_high[i]
    ))
  }
  cat(sprintf(
    "[OK]   Barrier topology model complete (%d/%d bootstrap reps successful). Wrote %s, %s\n",
    n_boot_success, N_BOOTSTRAP,
    file.path(output_dir_data, "barrier_topology_model_results.parquet"),
    file.path(output_dir_reports, "barrier_topology_models.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
