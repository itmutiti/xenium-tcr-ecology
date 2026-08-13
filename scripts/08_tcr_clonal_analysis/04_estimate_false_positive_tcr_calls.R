#!/usr/bin/env Rscript
# `08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R` -- 04_estimate_false_positive_tcr_calls.R
#
# Uses off-patient probes, non-T cells and spatial permutation as
# empirical negative controls for `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s cell-level TCR detection
# calls, restricted to `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s 105 identified probes.
#
# Three independent controls (see
# scripts/08_tcr_clonal_analysis/_04_prepare_false_positive_inputs.py and
# src/xenium_tcr_ecology/tcr/false_positive_estimation.py for how each is
# computed):
#   1. off_patient_detection_rate -- the same probe's detection rate in
#      the other candidate patients' own T cells (same manufacturing
#      batch, `08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`, `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`, not the intended patient).
#   2. non_tcell_detection_rate -- the probe's detection rate among its
#      own intended patient's non-T cells.
#   3. spatial_morans_i -- spatial autocorrelation of detection status
#      among the intended patient's own T cells (clonal signal should
#      cluster spatially; pure noise should not, the same logic already
#      validated in `07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py`).
#
# A per-probe empirical false-positive-rate estimate is computed as
# background_rate (the mean of the two rate-based controls) scaled by the
# probe's own number of evaluated T cells, expressed as a fraction of its
# own observed positive calls -- i.e. "what fraction of my own positive
# calls could plausibly be background-level noise alone."
#
# This script REPORTS the estimate; it does not itself exclude/re-call
# any `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py` detections -- that is `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`, `08_tcr_clonal_analysis/08_generate_tcr_release_report.py`'s job, using this
# script's output as one input among others.
#
# Primary output: reports/tcr/false_positive_model.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

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

compute_background_rate <- function(off_patient_rate, non_tcell_rate) {
  # Mean of the two rate-based controls -- distinct negative-control
  # populations (a different patient's T cells; this patient's own
  # non-T cells). These are correlated across probes, not independent:
  # Pearson r=0.528 (p=7.3e-9), Spearman rho=0.641 (p=1.7e-13), n=105
  # (checked against data/derived/tcr_false_positive_
  # controls.parquet). Consistent with both partly
  # reflecting a shared per-probe background/specificity level (a
  # less-specific probe shows elevated background in both populations
  # together), not two fully independent noise sources. The mean is
  # still a reasonable point-estimate summary of the two; no variance/CI
  # is computed on background_rate or empirical_fpr anywhere downstream,
  # so this correlation does not currently bias any reported number --
  # it would need to be accounted for if uncertainty quantification were
  # ever added here.
  (off_patient_rate + non_tcell_rate) / 2
}

compute_empirical_fpr <- function(own_n_tcells, own_n_detected, background_rate) {
  # Expected false positives under the background rate, expressed as a
  # fraction of the probe's own observed positive calls -- capped at 1.0
  # (a background rate that would predict more false positives than were
  # actually observed is reported as "up to 100% could be background,"
  # not an impossible rate above 100%).
  expected_false_positives <- own_n_tcells * background_rate
  fpr <- expected_false_positives / own_n_detected
  pmin(fpr, 1.0)
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  helper_script <- file.path(project_root, "scripts", "08_tcr_clonal_analysis", "_04_prepare_false_positive_inputs.py")
  cat("[INFO] Running Python false-positive-control preparation helper...\n")
  exit_code <- system2("python3", c(shQuote(helper_script), "--project-root", shQuote(project_root)))
  if (exit_code != 0) {
    stop("Python false-positive-control preparation helper failed -- see its output above.")
  }

  controls <- as.data.frame(read_parquet(file.path(project_root, "data", "derived", "tcr_false_positive_controls.parquet")))

  controls$background_rate <- compute_background_rate(controls$off_patient_detection_rate, controls$non_tcell_detection_rate)
  controls$empirical_fpr <- compute_empirical_fpr(controls$own_n_tcells, controls$own_n_detected, controls$background_rate)
  controls$signal_to_noise_ratio <- controls$own_detection_rate / controls$background_rate
  controls$spatially_corroborated <- !is.na(controls$spatial_pval_norm) & controls$spatial_pval_norm < 0.05

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(controls, file.path(output_dir_data, "tcr_false_positive_estimates.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "tcr")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Mixed 2x2 grid (A, B on top; C on bottom-left; bottom-right left
  # blank by compose_panels()'s row-major fill for n=3), not a full
  # vertical stack or full horizontal row: a word processor auto-scales
  # a pasted image to fit the page's text width (~6.5in), shrinking
  # every dimension including font size by that ratio, so an
  # all-horizontal 3-in-a-row image shrinks font size the most; an
  # all-vertical 3-stack becomes excessively long. The 2x2 compromise
  # keeps the width (and therefore the shrink ratio) moderate while
  # keeping each panel roughly as tall as it needs to be.
  open_publication_pdf(file.path(output_dir_reports, "false_positive_model.pdf"), width = 16.0, height = 13.5)

  p1 <- ggplot(controls, aes(x = background_rate, y = own_detection_rate)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    geom_point(aes(colour = spatially_corroborated), size = 2.8, alpha = 0.75) +
    scale_colour_manual(values = c("FALSE" = PUB_COLORS$not_significant, "TRUE" = PUB_COLORS$primary_analysis)) +
    scale_x_log10() + scale_y_log10() +
    labs(
      subtitle = "Own vs. background detection rate, per probe",
      x = "Background detection rate (log10)", y = "Own detection rate (log10)", colour = "Spatially\ncorroborated"
    ) +
    theme_publication()

  p2 <- ggplot(controls, aes(x = empirical_fpr)) +
    geom_histogram(bins = 30, fill = PUB_COLORS$primary_analysis, colour = "white", linewidth = 0.3) +
    labs(
      subtitle = "Per-probe empirical false-positive-rate estimate",
      x = "Empirical FPR estimate", y = "Number of probes"
    ) +
    theme_publication()

  p3 <- ggplot(controls[!is.na(controls$spatial_morans_i), ], aes(x = spatial_morans_i)) +
    geom_histogram(bins = 20, fill = PUB_COLORS$sensitivity_analysis, colour = "white", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    labs(
      subtitle = sprintf("Spatial autocorrelation of detection status (%d/%d probes)", sum(!is.na(controls$spatial_morans_i)), nrow(controls)),
      x = "Moran's I", y = "Number of probes"
    ) +
    theme_publication()

  compose_panels(list(p1, p2, p3), ncol = 2)
  dev.off()

  cat("[OK]   False-positive model complete.\n")
  cat(sprintf(
    "[OK]   %d probes: median empirical FPR %.4f, median signal-to-noise ratio %.2fx, %d/%d spatially corroborated.\n",
    nrow(controls), stats::median(controls$empirical_fpr), stats::median(controls$signal_to_noise_ratio),
    sum(controls$spatially_corroborated), sum(!is.na(controls$spatial_pval_norm))
  ))
  cat(sprintf(
    "[OK]   Wrote %s, %s\n",
    file.path(output_dir_reports, "false_positive_model.pdf"),
    file.path(output_dir_data, "tcr_false_positive_estimates.parquet")
  ))
}

if (sys.nframe() == 0) {
  main()
}
