#!/usr/bin/env Rscript
# `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R` -- 03_compare_cellular_composition_patient_level.R
#
# Compares per-lineage cellular composition between this project's
# validated HPV-positive (n=4) and HPV-negative (n=4) patient groups
# (`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`'s registered `hpv_validated_positive_vs_
# negative` contrast) using patient as the sampling unit -- not cells,
# and not sections either: a patient with 2 sections contributes one
# per-patient value (an equal-weight average across that patient's
# sections), so a 2-section patient cannot out-vote a 1-section patient.
# `15_hpv_stratified_analysis/02_run_prospective_power_simulation.R`'s finding (n=4 vs n=4 requires Cohen's d=2.5 for 80% power;
# exact rank-test minimum achievable p-value = 0.0286) is why this
# script uses an exact Wilcoxon rank-sum test, not a normal-theory
# approximation.
#
# Primary output: reports/hpv/composition_models.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(yaml)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"
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

parse_project_root_arg <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- which(args == "--project-root")
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  NULL
}

#' Per-cell `section_id`, assigned by matching each `cell_id` against
#' the one `valid_section_ids` entry that is its prefix (cell_id format
#' is `{section_id}_{barcode}`, and this project's section_id values
#' never prefix-conflict with each other). Stops with an error on zero
#' or multiple matches rather than guessing.
assign_section_id_from_cell_id <- function(cell_ids, valid_section_ids) {
  vapply(cell_ids, function(cid) {
    matches <- valid_section_ids[vapply(valid_section_ids, function(s) startsWith(cid, s), logical(1))]
    if (length(matches) != 1) stop(sprintf("cell_id '%s' matched %d section_id prefix(es), expected exactly 1.", cid, length(matches)))
    matches
  }, character(1), USE.NAMES = FALSE)
}

#' Per-(section, lineage) fraction of that section's cells, summing to
#' 1.0 per section across `all_lineages`.
compute_section_lineage_fractions <- function(lineage, section_id, all_lineages) {
  counts <- table(factor(section_id), factor(lineage, levels = all_lineages))
  totals <- rowSums(counts)
  fractions <- counts / totals
  result <- as.data.frame(as.table(fractions))
  colnames(result) <- c("section_id", "lineage", "fraction")
  result
}

#' Per-(patient, lineage) fraction -- an equal-weight average across
#' that patient's sections (a patient with 2 sections is not weighted
#' twice), matching this milestone's "patient as the sampling unit"
#' requirement.
aggregate_patient_lineage_fractions <- function(section_fractions, section_to_patient) {
  section_fractions$patient_id <- section_to_patient[as.character(section_fractions$section_id)]
  stats::aggregate(fraction ~ patient_id + lineage, data = section_fractions, FUN = mean)
}

#' Per-lineage exact two-sided Wilcoxon rank-sum comparison between the
#' positive and negative patient groups.
compare_hpv_groups_by_lineage <- function(patient_fractions, positive_ids, negative_ids) {
  lineages <- sort(unique(patient_fractions$lineage))
  rows <- list()
  for (lineage_name in lineages) {
    sub <- patient_fractions[patient_fractions$lineage == lineage_name, ]
    pos_values <- sub$fraction[sub$patient_id %in% positive_ids]
    neg_values <- sub$fraction[sub$patient_id %in% negative_ids]
    test <- stats::wilcox.test(pos_values, neg_values, exact = TRUE)
    rows[[lineage_name]] <- data.frame(
      lineage = lineage_name,
      n_positive = length(pos_values),
      n_negative = length(neg_values),
      median_positive = stats::median(pos_values),
      median_negative = stats::median(neg_values),
      pvalue = test$p.value
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result$pvalue_bh <- stats::p.adjust(result$pvalue, method = "BH")
  result
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  contrasts_path <- file.path(project_root, "governance", "hpv_primary_contrasts.yaml")
  annotations_path <- file.path(project_root, "data", "derived", "final_cell_annotations.parquet")
  sample_manifest_path <- file.path(project_root, "metadata", "sample_manifest.tsv")

  for (p in c(contrasts_path, annotations_path, sample_manifest_path)) {
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
  contrast_sections <- sample_manifest$section_id[sample_manifest$patient_id %in% contrast_patients]
  section_to_patient <- setNames(sample_manifest$patient_id, sample_manifest$section_id)[contrast_sections]

  annotations <- as.data.frame(read_parquet(annotations_path))
  colnames(annotations)[colnames(annotations) == "__index_level_0__"] <- "cell_id"
  annotations$section_id <- assign_section_id_from_cell_id(annotations$cell_id, unique(sample_manifest$section_id))
  contrast_cells <- annotations[annotations$section_id %in% contrast_sections, ]

  cat(sprintf("[INFO] Data: n=%d cells across %d sections for the contrast's %d patients.\n", nrow(contrast_cells), length(contrast_sections), length(contrast_patients)))

  all_lineages <- sort(unique(annotations$final_lineage))
  section_fractions <- compute_section_lineage_fractions(contrast_cells$final_lineage, contrast_cells$section_id, all_lineages)
  patient_fractions <- aggregate_patient_lineage_fractions(section_fractions, section_to_patient)

  result <- compare_hpv_groups_by_lineage(patient_fractions, positive_ids, negative_ids)

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "hpv_composition_comparison_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "hpv")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(output_dir_reports, "composition_models.pdf"), width = 10, height = 7)

  plot_data <- result
  plot_data$lineage <- factor(plot_data$lineage, levels = plot_data$lineage[order(plot_data$pvalue)])
  p1 <- ggplot(plot_data, aes(x = lineage)) +
    geom_point(aes(y = median_positive), color = "firebrick", size = 3) +
    geom_point(aes(y = median_negative), color = "steelblue", size = 3) +
    labs(
      title = "Q(HPV): per-lineage composition, HPV-positive (red) vs HPV-negative (blue)",
      subtitle = sprintf("Patient-level medians, n=%d vs n=%d; exact Wilcoxon p-values, BH-adjusted across %d lineages", length(positive_ids), length(negative_ids), nrow(result)),
      x = NULL, y = "Per-patient lineage fraction (median)"
    ) +
    coord_flip() +
    theme_minimal()
  print(p1)

  dev.off()

  cat("[INFO] Per-lineage comparison (patient-level, exact Wilcoxon):\n")
  for (i in order(result$pvalue)) {
    cat(sprintf("[INFO]   %-28s pos_median=%.3f neg_median=%.3f p=%.4f p_bh=%.4f\n", result$lineage[i], result$median_positive[i], result$median_negative[i], result$pvalue[i], result$pvalue_bh[i]))
  }
  n_significant_bh <- sum(result$pvalue_bh < 0.05)
  cat(sprintf(
    "[OK]   Composition comparison complete: %d/%d lineage(s) significant (BH q<0.05). Wrote %s, %s\n",
    n_significant_bh, nrow(result),
    file.path(output_dir_data, "hpv_composition_comparison_results.parquet"),
    file.path(output_dir_reports, "composition_models.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
