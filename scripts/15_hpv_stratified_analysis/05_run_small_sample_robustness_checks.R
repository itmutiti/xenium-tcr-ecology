#!/usr/bin/env Rscript
# `15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R` -- 05_run_small_sample_robustness_checks.R
#
# Fragility check on the 3 most notable results from `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`, `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R` (not all 25 tested categories --
# running LOPO/bootstrap/exhaustive-permutation on every category would
# be disproportionate for an "expose fragility" sanity check): `T_cell`
# and `Erythroid` lineage fractions (`15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R` -- `Erythroid` hit the
# exact-test floor raw p=0.0286; `T_cell` showed a biologically
# plausible but non-significant trend), and `clone_structure`'s
# per-patient `ecological_structure_score` (`15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R` -- the largest
# apparent median gap of anything tested, also non-significant).
#
# Three complementary fragility checks:
#   - Exhaustive permutation test on the raw per-patient values (not
#     ranks) -- an independent cross-check against Phase 15.03/15.04's
#     rank-based exact Wilcoxon results.
#   - Leave-one-patient-out (LOPO): recomputes the exact test with each
#     of the 8 contrast patients removed in turn, flagging whether any
#     single patient's removal flips direction or significance.
#   - Percentile bootstrap CI on the median difference (resampling with
#     replacement within each 4-patient group) -- coarse-resolution
#     given only 4 patients per group (only C(7,4)=35 distinct
#     within-group resample compositions are possible), not presented
#     as if it were a large-sample bootstrap.
#
# Primary output: reports/hpv/robustness.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(yaml)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"
CONTRAST_ID <- "hpv_validated_positive_vs_negative"
N_BOOTSTRAP <- 2000

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

#' Exhaustive two-sided permutation test on raw values (statistic =
#' |mean(positive) - mean(negative)|), enumerating every way to split
#' `values` into a positive-sized and negative-sized group.
run_exhaustive_permutation_test <- function(values, labels, positive_label, negative_label) {
  n_positive <- sum(labels == positive_label)
  n_total <- length(values)
  observed <- abs(mean(values[labels == positive_label]) - mean(values[labels == negative_label]))

  index_combinations <- utils::combn(n_total, n_positive)
  n_permutations <- ncol(index_combinations)
  null_stats <- apply(index_combinations, 2, function(idx) {
    abs(mean(values[idx]) - mean(values[-idx]))
  })

  list(
    observed = observed,
    n_permutations = n_permutations,
    pvalue = mean(null_stats >= observed - 1e-12)
  )
}

#' Leave-one-patient-out sensitivity -- drops each patient in turn and
#' recomputes the exact Wilcoxon test and median difference on the
#' remaining patients.
run_leave_one_out_sensitivity <- function(patient_values, positive_ids, negative_ids) {
  full_pos <- patient_values$value[patient_values$patient_id %in% positive_ids]
  full_neg <- patient_values$value[patient_values$patient_id %in% negative_ids]
  full_direction <- sign(stats::median(full_pos) - stats::median(full_neg))

  rows <- list()
  for (patient_id in patient_values$patient_id) {
    remaining <- patient_values[patient_values$patient_id != patient_id, ]
    pos_ids <- setdiff(positive_ids, patient_id)
    neg_ids <- setdiff(negative_ids, patient_id)
    pos_values <- remaining$value[remaining$patient_id %in% pos_ids]
    neg_values <- remaining$value[remaining$patient_id %in% neg_ids]
    test <- stats::wilcox.test(pos_values, neg_values, exact = TRUE)
    direction <- sign(stats::median(pos_values) - stats::median(neg_values))
    rows[[patient_id]] <- data.frame(
      patient_removed = patient_id,
      n_positive = length(pos_values),
      n_negative = length(neg_values),
      median_difference = stats::median(pos_values) - stats::median(neg_values),
      pvalue = test$p.value,
      same_direction_as_full = (direction == full_direction)
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Percentile bootstrap CI (resampling with replacement within each
#' group) for the median difference.
run_bootstrap_median_difference_ci <- function(patient_values, positive_ids, negative_ids, n_bootstrap = N_BOOTSTRAP) {
  pos_values <- patient_values$value[patient_values$patient_id %in% positive_ids]
  neg_values <- patient_values$value[patient_values$patient_id %in% negative_ids]
  point_estimate <- stats::median(pos_values) - stats::median(neg_values)

  boot_diffs <- vapply(seq_len(n_bootstrap), function(i) {
    pos_resample <- sample(pos_values, length(pos_values), replace = TRUE)
    neg_resample <- sample(neg_values, length(neg_values), replace = TRUE)
    stats::median(pos_resample) - stats::median(neg_resample)
  }, numeric(1))

  ci <- stats::quantile(boot_diffs, c(0.025, 0.975))
  list(point_estimate = point_estimate, ci_low = as.numeric(ci[1]), ci_high = as.numeric(ci[2]))
}

#' Recomputes the per-patient value vector for one of the 3 outcomes
#' this milestone stress-tests, reusing the aggregation logic
#' `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`, `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R` already established (redefined
#' here per this project's R-script-independence convention).
compute_patient_outcome_values <- function(outcome, project_root, section_to_patient) {
  if (outcome %in% c("T_cell", "Erythroid")) {
    annotations <- as.data.frame(read_parquet(file.path(project_root, "data", "derived", "final_cell_annotations.parquet")))
    colnames(annotations)[colnames(annotations) == "__index_level_0__"] <- "cell_id"
    valid_sections <- names(section_to_patient)
    annotations$section_id <- vapply(annotations$cell_id, function(cid) {
      matches <- valid_sections[vapply(valid_sections, function(s) startsWith(cid, s), logical(1))]
      matches[1]
    }, character(1))
    all_lineages <- sort(unique(annotations$final_lineage))
    counts <- table(factor(annotations$section_id), factor(annotations$final_lineage, levels = all_lineages))
    totals <- rowSums(counts)
    fractions <- as.data.frame(as.table(counts / totals))
    colnames(fractions) <- c("section_id", "lineage", "fraction")
    sub <- fractions[fractions$lineage == outcome, ]
    sub$patient_id <- section_to_patient[as.character(sub$section_id)]
    patient_means <- stats::aggregate(fraction ~ patient_id, data = sub, FUN = mean)
    return(data.frame(patient_id = patient_means$patient_id, value = patient_means$fraction))
  }
  if (outcome == "clone_structure") {
    structure <- as.data.frame(read_parquet(file.path(project_root, "data", "releases", "v1_clone_structure", "clone_ecological_structure.parquet")))
    section_means <- stats::aggregate(ecological_structure_score ~ section_id, data = structure, FUN = mean)
    section_means$patient_id <- section_to_patient[as.character(section_means$section_id)]
    patient_means <- stats::aggregate(ecological_structure_score ~ patient_id, data = section_means, FUN = mean)
    return(data.frame(patient_id = patient_means$patient_id, value = patient_means$ecological_structure_score))
  }
  stop(sprintf("Unknown outcome '%s'.", outcome))
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  contrasts_path <- file.path(project_root, "governance", "hpv_primary_contrasts.yaml")
  sample_manifest_path <- file.path(project_root, "metadata", "sample_manifest.tsv")
  for (p in c(contrasts_path, sample_manifest_path)) {
    if (!file.exists(p)) stop(sprintf("'%s' not found.", p))
  }

  contrasts <- yaml::read_yaml(contrasts_path)
  contrast <- Filter(function(c) c$contrast_id == CONTRAST_ID, contrasts$primary_contrasts)
  if (length(contrast) == 0) stop(sprintf("Contrast '%s' not found in '%s'. Run `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py` first.", CONTRAST_ID, contrasts_path))
  contrast <- contrast[[1]]
  positive_ids <- unlist(contrast$positive_group$patient_ids)
  negative_ids <- unlist(contrast$negative_group$patient_ids)

  sample_manifest <- read.delim(sample_manifest_path, sep = "\t", stringsAsFactors = FALSE)
  section_to_patient <- setNames(sample_manifest$patient_id, sample_manifest$section_id)

  outcomes <- c("T_cell", "Erythroid", "clone_structure")
  lopo_rows <- list()
  summary_rows <- list()

  set.seed(RNG_SEED)
  for (outcome in outcomes) {
    patient_values <- compute_patient_outcome_values(outcome, project_root, section_to_patient)
    labels <- ifelse(patient_values$patient_id %in% positive_ids, "positive", ifelse(patient_values$patient_id %in% negative_ids, "negative", NA))
    patient_values <- patient_values[!is.na(labels), ]
    labels <- labels[!is.na(labels)]

    perm_result <- run_exhaustive_permutation_test(patient_values$value, labels, "positive", "negative")
    lopo_result <- run_leave_one_out_sensitivity(patient_values, positive_ids, negative_ids)
    lopo_result$outcome <- outcome
    lopo_rows[[outcome]] <- lopo_result
    boot_result <- run_bootstrap_median_difference_ci(patient_values, positive_ids, negative_ids)

    summary_rows[[outcome]] <- data.frame(
      outcome = outcome,
      n_permutations = perm_result$n_permutations,
      permutation_pvalue = perm_result$pvalue,
      bootstrap_point_estimate = boot_result$point_estimate,
      bootstrap_ci_low = boot_result$ci_low,
      bootstrap_ci_high = boot_result$ci_high,
      n_lopo_direction_flips = sum(!lopo_result$same_direction_as_full)
    )
  }

  summary_result <- do.call(rbind, summary_rows)
  rownames(summary_result) <- NULL
  lopo_all <- do.call(rbind, lopo_rows)
  rownames(lopo_all) <- NULL

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(summary_result, file.path(output_dir_data, "hpv_robustness_summary.parquet"))
  write_parquet(lopo_all, file.path(output_dir_data, "hpv_robustness_lopo.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "hpv")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # A on top, B (a 3-outcome facet grid) below; B's facets run in a
  # single row (ncol=3) rather than stacked, keeping the combined
  # canvas close to square instead of excessively tall.
  open_publication_pdf(file.path(output_dir_reports, "robustness.pdf"), width = 16.5, height = 10.0)

  p1 <- ggplot(summary_result, aes(x = outcome, y = bootstrap_point_estimate)) +
    geom_col(fill = PUB_COLORS$sensitivity_analysis, width = 0.6) +
    geom_errorbar(aes(ymin = bootstrap_ci_low, ymax = bootstrap_ci_high), width = 0.25, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    labs(
      subtitle = sprintf("Bootstrap CI on median difference (n = 4/group, %d replicates)", N_BOOTSTRAP),
      x = NULL, y = "Median difference\n(positive - negative)"
    ) +
    theme_publication()

  p2 <- ggplot(lopo_all, aes(x = patient_removed, y = median_difference, color = same_direction_as_full)) +
    geom_point(size = 3.4) +
    facet_wrap(~outcome, scales = "free_y", ncol = 3) +
    scale_color_manual(values = setNames(c(PUB_COLORS$primary_analysis, PUB_COLORS$not_significant), c(TRUE, FALSE))) +
    labs(
      subtitle = "Median difference after removing each patient in turn",
      x = "Patient removed", y = "Median difference", color = "Same direction\nas full sample"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing.x = unit(1.4, "lines"))

  compose_panels(list(p1, p2), ncol = 1, heights = c(1, 1.3))
  dev.off()

  cat("[INFO] Robustness summary:\n")
  for (i in seq_len(nrow(summary_result))) {
    cat(sprintf(
      "[INFO]   %-16s perm_p=%.4f (n_perm=%d) bootstrap=%.4f (%.4f-%.4f) lopo_direction_flips=%d/8\n",
      summary_result$outcome[i], summary_result$permutation_pvalue[i], summary_result$n_permutations[i],
      summary_result$bootstrap_point_estimate[i], summary_result$bootstrap_ci_low[i], summary_result$bootstrap_ci_high[i],
      summary_result$n_lopo_direction_flips[i]
    ))
  }
  cat(sprintf(
    "[OK]   Robustness checks complete. Wrote %s, %s, %s\n",
    file.path(output_dir_data, "hpv_robustness_summary.parquet"),
    file.path(output_dir_data, "hpv_robustness_lopo.parquet"),
    file.path(output_dir_reports, "robustness.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
