#!/usr/bin/env Rscript
# `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R` -- 04_compare_ecosystem_and_clone_structure.R
#
# Tests `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`'s registered `hpv_validated_positive_vs_
# negative` contrast (n=4 vs n=4) against two outcome domains: (1)
# ecosystem abundance/topology (`10_niche_and_ecosystem_discovery/05_quantify_ecosystem_abundance_and_topology.py`'s `ecosystem_metrics.
# parquet` -- 6 ecosystem labels x `abundance` + `mixing_index`), and
# (2) this project's frozen, confirmed clone ecological-structure score
# (`11_clone_spatial_descriptors/06_discover_provisional_structure.R`, Clone Ecology Confirmatory Models's `clone_ecological_structure.
# parquet`, taxonomy_version=v1_provisional). Patient-as-sampling-unit
# discipline and exact Wilcoxon tests, matching Phase 15.03's
# convention (this script does not source 15.03 -- the small
# comparison/aggregation logic is independently redefined here, per
# this project's convention of not cross-sourcing between R scripts). A
# single BH correction is applied across all 13 tests run by this one
# script (6 ecosystem labels x 2 metrics + 1 clone-structure test),
# since they are one combined test family, not 13 independent analyses.
#
# Primary output: reports/hpv/structure_models.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(yaml)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"
CONTRAST_ID <- "hpv_validated_positive_vs_negative"
ECOSYSTEM_METRICS <- c("abundance", "mixing_index")

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

#' Per-(patient, category) value -- an equal-weight average across
#' that patient's sections (matching Phase 15.03's "patient as the
#' sampling unit" convention).
aggregate_patient_category_metric <- function(section_id, category, value, section_to_patient) {
  patient_id <- section_to_patient[as.character(section_id)]
  stats::aggregate(value ~ patient_id + category, data = data.frame(patient_id = patient_id, category = category, value = value), FUN = mean)
}

#' Per-patient clone ecological-structure score, collapsed in two
#' equal-weight stages -- first the mean across a section's clones,
#' then the mean across a patient's sections -- so a section with many
#' clones cannot out-vote a section with few, and a patient with many
#' sections cannot out-vote a patient with few.
aggregate_patient_clone_structure <- function(section_id, ecological_structure_score, section_to_patient) {
  section_means <- stats::aggregate(ecological_structure_score ~ section_id, data = data.frame(section_id = section_id, ecological_structure_score = ecological_structure_score), FUN = mean)
  section_means$patient_id <- section_to_patient[as.character(section_means$section_id)]
  patient_means <- stats::aggregate(ecological_structure_score ~ patient_id, data = section_means, FUN = mean)
  data.frame(patient_id = patient_means$patient_id, value = patient_means$ecological_structure_score)
}

#' Exact two-sided Wilcoxon rank-sum comparison between positive and
#' negative patient groups for a single `category` (or the whole
#' `patient_values`, if it has no `category` column with >1 level).
compare_hpv_groups <- function(patient_values, positive_ids, negative_ids) {
  pos_values <- patient_values$value[patient_values$patient_id %in% positive_ids]
  neg_values <- patient_values$value[patient_values$patient_id %in% negative_ids]
  test <- stats::wilcox.test(pos_values, neg_values, exact = TRUE)
  data.frame(
    n_positive = length(pos_values),
    n_negative = length(neg_values),
    median_positive = stats::median(pos_values),
    median_negative = stats::median(neg_values),
    pvalue = test$p.value
  )
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  contrasts_path <- file.path(project_root, "governance", "hpv_primary_contrasts.yaml")
  ecosystem_path <- file.path(project_root, "data", "derived", "ecosystem_metrics.parquet")
  structure_path <- file.path(project_root, "data", "releases", "v1_clone_structure", "clone_ecological_structure.parquet")
  sample_manifest_path <- file.path(project_root, "metadata", "sample_manifest.tsv")

  for (p in c(contrasts_path, ecosystem_path, structure_path, sample_manifest_path)) {
    if (!file.exists(p)) stop(sprintf("'%s' not found.", p))
  }

  contrasts <- yaml::read_yaml(contrasts_path)
  contrast <- Filter(function(c) c$contrast_id == CONTRAST_ID, contrasts$primary_contrasts)
  if (length(contrast) == 0) stop(sprintf("Contrast '%s' not found in '%s'. Run `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py` first.", CONTRAST_ID, contrasts_path))
  contrast <- contrast[[1]]
  positive_ids <- unlist(contrast$positive_group$patient_ids)
  negative_ids <- unlist(contrast$negative_group$patient_ids)

  cat(sprintf("[INFO] Contrast '%s': n=%d positive vs n=%d negative.\n", CONTRAST_ID, length(positive_ids), length(negative_ids)))

  sample_manifest <- read.delim(sample_manifest_path, sep = "\t", stringsAsFactors = FALSE)
  contrast_patients <- c(positive_ids, negative_ids)
  section_to_patient <- setNames(sample_manifest$patient_id, sample_manifest$section_id)

  ecosystem <- as.data.frame(read_parquet(ecosystem_path))
  ecosystem <- ecosystem[section_to_patient[ecosystem$section_id] %in% contrast_patients, ]

  structure <- as.data.frame(read_parquet(structure_path))
  structure <- structure[structure$patient_id %in% contrast_patients, ]
  if (nrow(structure) == 0) stop("No clone-structure rows found for the contrast's patients.")

  rows <- list()
  for (metric in ECOSYSTEM_METRICS) {
    patient_values <- aggregate_patient_category_metric(ecosystem$section_id, ecosystem$ecosystem_label, ecosystem[[metric]], section_to_patient)
    for (label in sort(unique(patient_values$category))) {
      sub <- patient_values[patient_values$category == label, c("patient_id", "value")]
      comparison <- compare_hpv_groups(sub, positive_ids, negative_ids)
      rows[[paste(metric, label)]] <- cbind(outcome_domain = "ecosystem", category = label, metric = metric, comparison)
    }
  }

  clone_structure_patient_values <- aggregate_patient_clone_structure(structure$section_id, structure$ecological_structure_score, section_to_patient)
  structure_comparison <- compare_hpv_groups(clone_structure_patient_values, positive_ids, negative_ids)
  rows[["clone_structure"]] <- cbind(outcome_domain = "clone_structure", category = "ecological_structure_score", metric = "ecological_structure_score", structure_comparison)

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result$pvalue_bh <- stats::p.adjust(result$pvalue, method = "BH")

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "hpv_structure_comparison_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "hpv")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(output_dir_reports, "structure_models.pdf"), width = 10, height = 7)

  plot_data <- result
  plot_data$label <- paste(plot_data$outcome_domain, plot_data$category, plot_data$metric, sep = " | ")
  plot_data$label <- factor(plot_data$label, levels = plot_data$label[order(plot_data$pvalue)])
  p1 <- ggplot(plot_data, aes(x = label)) +
    geom_point(aes(y = median_positive), color = "firebrick", size = 3) +
    geom_point(aes(y = median_negative), color = "steelblue", size = 3) +
    labs(
      title = "Q(HPV): ecosystem/clone-structure outcomes, HPV-positive (red) vs HPV-negative (blue)",
      subtitle = sprintf("Patient-level medians, n=%d vs n=%d; exact Wilcoxon p-values, BH-adjusted across %d tests", length(positive_ids), length(negative_ids), nrow(result)),
      x = NULL, y = "Per-patient value (median)"
    ) +
    coord_flip() +
    theme_minimal()
  print(p1)

  dev.off()

  cat("[INFO] Ecosystem/clone-structure comparison (patient-level, exact Wilcoxon):\n")
  for (i in order(result$pvalue)) {
    cat(sprintf("[INFO]   %-12s %-32s %-24s pos=%.3f neg=%.3f p=%.4f p_bh=%.4f\n", result$outcome_domain[i], result$category[i], result$metric[i], result$median_positive[i], result$median_negative[i], result$pvalue[i], result$pvalue_bh[i]))
  }
  n_significant_bh <- sum(result$pvalue_bh < 0.05)
  cat(sprintf(
    "[OK]   Ecosystem/clone-structure comparison complete: %d/%d test(s) significant (BH q<0.05). Wrote %s, %s\n",
    n_significant_bh, nrow(result),
    file.path(output_dir_data, "hpv_structure_comparison_results.parquet"),
    file.path(output_dir_reports, "structure_models.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
