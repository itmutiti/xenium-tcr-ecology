#!/usr/bin/env Rscript
# `15_hpv_stratified_analysis/02_run_prospective_power_simulation.R` -- 02_run_prospective_power_simulation.R
#
# Simulates the minimum detectable effect size (MDES) for this project's
# validated primary HPV contrast -- n=4 HPV-positive vs n=4 HPV-negative
# (`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`'s `governance/hpv_primary_contrasts.yaml`, not the
# scaffold's originally-anticipated "n=5 vs n=5" -- `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`'s discordance
# finding reduced the usable n before this milestone ran; loaded
# programmatically from the registered contrast, not hardcoded).
#
# "Planned model" simulated: a patient-level, two-sample comparison of a
# continuous summary statistic (matching `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`, `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R`'s
# patient-as-sampling-unit design) -- both a two-sample t-test (the
# natural parametric baseline) and the exact minimum achievable
# two-sided p-value for a rank-based permutation/exact test at this n
# (`15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R` will use permutation/exact tests specifically
# because of the same small-n concern this milestone quantifies here).
#
# Primary output: reports/hpv/power_simulation.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

ALPHA <- 0.05
TARGET_POWER <- 0.8
N_SIMULATIONS <- 4000
EFFECT_SIZE_GRID <- seq(0.5, 5.0, by = 0.25)
CONTRAST_ID <- "hpv_validated_positive_vs_negative"

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

#' Pure, testable: Monte Carlo power for a two-sample t-test, given
#' group sizes and a standardised effect size (Cohen's d, equal-variance
#' unit-sd groups).
simulate_power_two_sample_t <- function(n1, n2, cohens_d, alpha, n_simulations) {
  rejections <- 0
  for (i in seq_len(n_simulations)) {
    group1 <- stats::rnorm(n1, mean = 0, sd = 1)
    group2 <- stats::rnorm(n2, mean = cohens_d, sd = 1)
    pvalue <- stats::t.test(group1, group2)$p.value
    if (pvalue < alpha) rejections <- rejections + 1
  }
  rejections / n_simulations
}

#' Pure, testable: exact (closed-form, not simulated) minimum achievable
#' two-sided p-value for a rank-based exact/permutation test at group
#' sizes `n1`, `n2` -- the smallest possible p-value is
#' `2 / choose(n1+n2, n1)`, achieved only when the two groups are
#' perfectly separated (every value in one group exceeds every value in
#' the other).
compute_exact_permutation_min_pvalue <- function(n1, n2) {
  2 / choose(n1 + n2, n1)
}

#' Pure, testable: smallest `effect_sizes` grid value whose simulated
#' power (`simulate_power_two_sample_t`) reaches `target_power`; `NA`
#' (not a silent extrapolation) if no grid value clears it.
find_minimum_detectable_effect <- function(n1, n2, alpha, target_power, effect_sizes, n_simulations) {
  powers <- vapply(effect_sizes, function(d) simulate_power_two_sample_t(n1, n2, d, alpha, n_simulations), numeric(1))
  cleared <- effect_sizes[powers >= target_power]
  list(
    minimum_detectable_effect = if (length(cleared) > 0) min(cleared) else NA_real_,
    effect_sizes = effect_sizes,
    powers = powers
  )
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  contrasts_path <- file.path(project_root, "governance", "hpv_primary_contrasts.yaml")
  if (!file.exists(contrasts_path)) stop(sprintf("'%s' not found. Run `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py` first.", contrasts_path))

  contrasts <- yaml::read_yaml(contrasts_path)
  contrast <- Filter(function(c) c$contrast_id == CONTRAST_ID, contrasts$primary_contrasts)
  if (length(contrast) == 0) stop(sprintf("Contrast '%s' not found in '%s'. Run `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py` first.", CONTRAST_ID, contrasts_path))
  contrast <- contrast[[1]]
  n1 <- contrast$n_positive
  n2 <- contrast$n_negative

  cat(sprintf("[INFO] Registered contrast '%s': n1=%d (positive) vs n2=%d (negative).\n", CONTRAST_ID, n1, n2))

  set.seed(RNG_SEED)
  mde_result <- find_minimum_detectable_effect(n1, n2, ALPHA, TARGET_POWER, EFFECT_SIZE_GRID, N_SIMULATIONS)
  exact_min_pvalue <- compute_exact_permutation_min_pvalue(n1, n2)

  power_curve <- data.frame(cohens_d = mde_result$effect_sizes, power = mde_result$powers)

  output_dir_data <- file.path(project_root, "data", "derived")
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(power_curve, file.path(output_dir_data, "hpv_power_simulation_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "hpv")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Narrower, taller canvas: a word processor auto-scales a pasted
  # image to fit the page's text width (~6.5in), shrinking every
  # dimension including font size by that ratio. A narrower canvas
  # shrinks less, so the same in-image point size reads larger once
  # pasted. Title/subtitle are wrapped onto separate lines below to
  # remain fully legible at this narrower width.
  open_publication_pdf(file.path(output_dir_reports, "power_simulation.pdf"), width = 9.6, height = 7.6)

  mde_label <- if (is.na(mde_result$minimum_detectable_effect)) {
    sprintf("No effect size up to d = %.2f reached %.0f%% power", max(EFFECT_SIZE_GRID), TARGET_POWER * 100)
  } else {
    sprintf("Minimum detectable effect (%.0f%% power): d = %.2f", TARGET_POWER * 100, mde_result$minimum_detectable_effect)
  }
  p1 <- ggplot(power_curve, aes(x = cohens_d, y = power)) +
    geom_line(color = PUB_COLORS$sensitivity_analysis, linewidth = 1.3) +
    geom_point(color = PUB_COLORS$sensitivity_analysis, size = 2.6) +
    geom_hline(yintercept = TARGET_POWER, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    labs(
      title = sprintf("Prospective power, n = %d vs n = %d\n(HPV primary contrast)", n1, n2),
      subtitle = sprintf("%s\nExact minimum achievable two-sided p-value at this n = %.4f", mde_label, exact_min_pvalue),
      x = "Cohen's d (standardised effect size)", y = sprintf("Simulated power\n(n = %d Monte Carlo reps/point)", N_SIMULATIONS)
    ) +
    ylim(0, 1) +
    theme_publication()
  print(p1)

  dev.off()

  # Backfills the computed MDES into the same registered contrast entry
  # in governance/hpv_primary_contrasts.yaml (left as a `null`
  # placeholder by `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`, whose module docstring states
  # this is `15_hpv_stratified_analysis/02_run_prospective_power_simulation.R`'s job) -- every other field of every other
  # contrast entry is preserved unchanged.
  for (i in seq_along(contrasts$primary_contrasts)) {
    if (contrasts$primary_contrasts[[i]]$contrast_id == CONTRAST_ID) {
      # Deliberate: assigning NA_real_ (not NULL) keeps the key present
      # in the written YAML even when no grid value cleared the target
      # power -- `x$field <- NULL` would silently delete the list
      # element in R, which would make a downstream reader unable to
      # distinguish "not yet computed" from "computed, found nothing".
      contrasts$primary_contrasts[[i]]$minimum_detectable_effect <- mde_result$minimum_detectable_effect
    }
  }
  writeLines(
    c(
      "# Written by scripts/15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py (`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`), before any HPV",
      "# hypothesis test is run. A maximum of two entries are permitted here for this project;",
      "# every other HPV comparison anywhere in the pipeline is exploratory only and must not reference",
      "# hpv_status in a confirmatory model formula.",
      "# minimum_detectable_effect is backfilled by scripts/15_hpv_stratified_analysis/02_run_prospective_power_simulation.R.",
      yaml::as.yaml(contrasts)
    ),
    contrasts_path
  )

  cat(sprintf("[INFO] Exact rank-test minimum achievable two-sided p-value at n=%d vs n=%d: %.4f\n", n1, n2, exact_min_pvalue))
  cat(sprintf("[INFO] %s\n", mde_label))
  cat(sprintf(
    "[OK]   Power simulation complete. Wrote %s, %s\n",
    file.path(output_dir_data, "hpv_power_simulation_results.parquet"),
    file.path(output_dir_reports, "power_simulation.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
