#!/usr/bin/env Rscript
# `04_quality_control/06_define_qc_thresholds_hierarchically.R` -- 06_define_qc_thresholds_hierarchically.R
#
# Defines hierarchical, section-aware cell-level QC thresholds: an absolute
# floor on transcript counts and genes detected, plus a per-section robust
# (median/MAD) outlier check on transcript counts (same modified z-score
# convention as `04_quality_control/02_detect_spatial_qc_artifacts.py`, threshold 3.5) so that sections with
# lower sequencing/imaging depth are not unfairly decimated by a single
# global cutoff. Also defines a ceiling on negative-control-probe/codeword
# ratio, and documents sensitivity alternatives (lenient/standard/strict)
# with exclusion-rate estimates against the actual QC metrics, rather
# than committing to a single unexamined cut-off.
#
# This script only *defines* thresholds; it does not exclude any cells or
# transcripts -- `04_quality_control/07_apply_qc_filters_with_audit_trail.py` applies them, with a reason code recorded for
# every excluded cell (config/qc_thresholds.yaml is 4.07's input).
#
# The "standard" profile was reviewed against per-section and pooled
# exclusion-rate estimates and selected on 2026-07-10.
#
# Primary output: config/qc_thresholds.yaml

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
  library(yaml)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

# Mirrors the precedence and marker-file convention of the Python
# implementation (src/xenium_tcr_ecology/infra/paths.py::find_project_root):
# --project-root CLI arg, then the environment variable, then walking up
# from this script's location looking for the marker file.
find_project_root <- function(cli_arg = NULL) {
  check_marker <- function(candidate) file.exists(file.path(candidate, PROJECT_ROOT_MARKER))

  if (!is.null(cli_arg) && nzchar(cli_arg)) {
    candidate <- normalizePath(cli_arg, mustWork = FALSE)
    if (check_marker(candidate)) return(candidate)
    stop(sprintf(
      "Project root supplied via --project-root ('%s') does not contain '%s'.",
      candidate, PROJECT_ROOT_MARKER
    ))
  }

  env_root <- Sys.getenv(ENV_ROOT_VAR, unset = NA)
  if (!is.na(env_root) && nzchar(env_root)) {
    candidate <- normalizePath(env_root, mustWork = FALSE)
    if (check_marker(candidate)) return(candidate)
    stop(sprintf(
      "Project root supplied via $%s ('%s') does not contain '%s'.",
      ENV_ROOT_VAR, candidate, PROJECT_ROOT_MARKER
    ))
  }

  here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
  parent <- here
  repeat {
    if (check_marker(parent)) return(parent)
    next_parent <- dirname(parent)
    if (next_parent == parent) break
    parent <- next_parent
  }

  stop(sprintf(
    "Could not locate the project root. None of --project-root, $%s, or a '%s' marker file found by walking up from this script identified it.",
    ENV_ROOT_VAR, PROJECT_ROOT_MARKER
  ))
}

parse_project_root_arg <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- which(args == "--project-root")
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  NULL
}

# Robust (median/MAD) modified z-score -- identical convention to
# `04_quality_control/02_detect_spatial_qc_artifacts.py`'s `_modified_z_scores()` (Iglewicz & Hoaglin, 1993).
modified_z_scores <- function(x) {
  med <- median(x)
  mad_val <- median(abs(x - med))
  if (mad_val == 0) return(rep(0, length(x)))
  0.6745 * (x - med) / mad_val
}

# Threshold profiles reviewed against exclusion-rate estimates and
# signed off on 2026-07-10.
# "standard" is the active default; "lenient" and "strict" are the
# documented sensitivity alternatives.
PROFILES <- list(
  lenient = list(
    min_transcript_counts = 5,
    min_genes_detected = 3,
    max_control_probe_ratio = 0.10,
    max_control_codeword_ratio = 0.10,
    section_relative_min_counts_z = -5.0
  ),
  standard = list(
    min_transcript_counts = 10,
    min_genes_detected = 3,
    max_control_probe_ratio = 0.05,
    max_control_codeword_ratio = 0.05,
    section_relative_min_counts_z = -3.5
  ),
  strict = list(
    min_transcript_counts = 20,
    min_genes_detected = 5,
    max_control_probe_ratio = 0.02,
    max_control_codeword_ratio = 0.02,
    section_relative_min_counts_z = -3.0
  )
)
ACTIVE_PROFILE <- "standard"

evaluate_profile <- function(df, profile) {
  df$transcript_counts < profile$min_transcript_counts |
    df$n_genes_detected < profile$min_genes_detected |
    df$control_probe_ratio > profile$max_control_probe_ratio |
    df$control_codeword_ratio > profile$max_control_codeword_ratio |
    df$z_counts < profile$section_relative_min_counts_z
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  cell_qc_path <- file.path(project_root, "data", "derived", "cell_qc_metrics.parquet")
  fov_path <- file.path(project_root, "reports", "qc", "spatial_artifact_masks", "fov_artifact_candidates.tsv")
  output_path <- file.path(project_root, "config", "qc_thresholds.yaml")

  if (!file.exists(cell_qc_path)) {
    stop(sprintf("'%s' not found. Run `04_quality_control/00_compute_cell_level_qc_metrics.py` first.", cell_qc_path))
  }
  if (!file.exists(fov_path)) {
    stop(sprintf("'%s' not found. Run `04_quality_control/02_detect_spatial_qc_artifacts.py` first.", fov_path))
  }

  df <- read_parquet(cell_qc_path)
  fov_df <- read_tsv(fov_path, show_col_types = FALSE)

  df <- df %>%
    group_by(section_id) %>%
    mutate(z_counts = modified_z_scores(transcript_counts)) %>%
    ungroup()

  n_total <- nrow(df)

  # Sensitivity analysis: pooled and per-section exclusion counts under each
  # profile, computed against the data -- not asserted, measured.
  pooled_sensitivity <- list()
  per_section_sensitivity <- list()
  for (name in names(PROFILES)) {
    excluded <- evaluate_profile(df, PROFILES[[name]])
    pooled_sensitivity[[name]] <- list(
      n_excluded = as.integer(sum(excluded)),
      fraction_excluded = round(sum(excluded) / n_total, 4)
    )
    section_frac <- tapply(excluded, df$section_id, mean)
    per_section_sensitivity[[name]] <- as.list(round(section_frac, 4))
  }

  # No-nucleus and FOV-artifact flags are informational at this stage
  # (defined here, applied -- or not -- in `04_quality_control/07_apply_qc_filters_with_audit_trail.py`). A missing nucleus
  # segmentation does not by itself indicate the transcript signal for that
  # cell is unreliable, so this profile flags rather than excludes it.
  no_nucleus_fraction <- round(mean(is.na(df$nucleus_area)), 4)

  n_fovs_total <- nrow(fov_df)
  n_fovs_flagged <- sum(fov_df$flagged_artifact_candidate)
  transcript_fraction_in_flagged_fovs <- round(
    sum(fov_df$n_transcripts[fov_df$flagged_artifact_candidate]) / sum(fov_df$n_transcripts), 4
  )

  config <- list(
    active_profile = ACTIVE_PROFILE,
    generated_by = "scripts/04_quality_control/06_define_qc_thresholds_hierarchically.R",
    generated_date = as.character(Sys.Date()),
    rationale = paste(
      "Hierarchical / section-aware framework, not a single global cut-off:",
      "(1) an absolute floor on transcript_counts and n_genes_detected;",
      "(2) a per-section robust (median/MAD) outlier check on transcript_counts",
      "(same modified z-score convention as `04_quality_control/02_detect_spatial_qc_artifacts.py`, threshold configurable per profile),",
      "since sequencing/imaging depth varies meaningfully by section/batch;",
      "(3) a ceiling on negative-control-probe/codeword ratio;",
      "(4) FOV-level artifact exclusion applied by joining against `04_quality_control/02_detect_spatial_qc_artifacts.py`'s",
      "fov_artifact_candidates.tsv at filter-application time (`04_quality_control/07_apply_qc_filters_with_audit_trail.py`), not duplicated here.",
      "The 'standard' profile was reviewed against the exclusion-rate estimates below and",
      "selected on 2026-07-10."
    ),
    profiles = PROFILES,
    flags = list(
      no_nucleus_detected = list(
        definition = "nucleus_area is missing (no segmented nucleus for this cell)",
        action = "flag_only",
        pooled_fraction = no_nucleus_fraction
      ),
      fov_artifact_candidate = list(
        definition = "cell's FOV is flagged in reports/qc/spatial_artifact_masks/fov_artifact_candidates.tsv",
        action = "exclude",
        applied_at = "04_quality_control/07_apply_qc_filters_with_audit_trail.py",
        n_fovs_total = as.integer(n_fovs_total),
        n_fovs_flagged = as.integer(n_fovs_flagged),
        estimated_transcript_volume_fraction_affected = transcript_fraction_in_flagged_fovs
      )
    ),
    sensitivity_analysis = list(
      n_cells_total = as.integer(n_total),
      pooled = pooled_sensitivity,
      per_section = per_section_sensitivity
    )
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_yaml(config, output_path)

  cat(sprintf(
    "[OK]   Active profile '%s': %d / %d cells (%.2f%%) would be excluded pooled. Wrote %s\n",
    ACTIVE_PROFILE,
    pooled_sensitivity[[ACTIVE_PROFILE]]$n_excluded,
    n_total,
    100 * pooled_sensitivity[[ACTIVE_PROFILE]]$fraction_excluded,
    output_path
  ))
}

# Guards against execution when this file is `source()`d (e.g. by a test
# suite) rather than run directly via `Rscript`. sys.nframe() is 0 at top
# level for direct execution and >0 when sourced -- verified empirically,
# mirrors Python's `if __name__ == "__main__":` idiom.
if (sys.nframe() == 0) {
  main()
}
