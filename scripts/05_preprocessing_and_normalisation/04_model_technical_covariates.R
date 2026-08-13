#!/usr/bin/env Rscript
# `05_preprocessing_and_normalisation/04_model_technical_covariates.R` -- 04_model_technical_covariates.R
#
# Estimates the contribution of technical factors (sequencing/imaging depth,
# negative-control-probe burden) and section/run (nested within patient) to
# each of `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s 8 program scores, using a linear mixed model per
# program: score ~ depth_z + control_burden_z + (1|patient_id) + (1|section_id).
# section_id values are globally unique (e.g. "P09_run1"), so `(1|section_id)`
# already correctly represents "run nested within patient" without an
# explicit interaction term.
#
# This is a DIAGNOSTIC decomposition, not a correction: nothing is
# regressed out of the analysis matrix here, and patient variance is
# reported, not erased -- exactly the blueprint's own framing
# ("without erasing patient biology").
#
# Reads program_scores.parquet and cell_qc_metrics.parquet directly (both
# already plain parquet, unlike `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`'s analysis_ready.h5ad) -- no
# Python subprocess bridge needed for this script.
#
# Primary output: reports/preprocess/variance_partition.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(lme4)
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

PROGRAMS <- c(
  "cytotoxicity_score", "exhaustion_score", "activation_score", "interferon_score",
  "proliferation_score", "stress_score", "emt_score", "antigen_presentation_score"
)
SUBSAMPLE_PER_SECTION <- 500

# Manual variance decomposition (Nakagawa & Schielzeth 2013-style marginal/
# conditional variance split), computed directly from lme4's own outputs.
# Neither MuMIn nor variancePartition is part of this project's locked R
# environment (environment/conda/main.yml's deliberate incremental-locking
# convention) -- this is a small, from-scratch equivalent, not a
# reimplementation of lme4 itself.
decompose_variance <- function(model) {
  fixed_pred <- as.numeric(model.matrix(model) %*% lme4::fixef(model))
  var_fixed <- var(fixed_pred)

  vc <- as.data.frame(lme4::VarCorr(model))
  get_vc <- function(group) {
    v <- vc$vcov[vc$grp == group]
    if (length(v) == 0) 0 else v
  }
  var_patient <- get_vc("patient_id")
  var_section <- get_vc("section_id")
  var_residual <- get_vc("Residual")

  total <- var_fixed + var_patient + var_section + var_residual
  data.frame(
    component = c("fixed (depth + control burden)", "patient", "section (run nested in patient)", "residual"),
    variance = c(var_fixed, var_patient, var_section, var_residual),
    fraction = c(var_fixed, var_patient, var_section, var_residual) / total
  )
}

fit_and_decompose <- function(df, outcome) {
  formula_str <- sprintf("%s ~ depth_z + control_burden_z + (1|patient_id) + (1|section_id)", outcome)
  model <- lme4::lmer(as.formula(formula_str), data = df, REML = TRUE)
  result <- decompose_variance(model)
  result$program <- outcome
  result$converged <- length(model@optinfo$conv$lme4$messages) == 0
  result
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  scores_path <- file.path(project_root, "data", "derived", "program_scores.parquet")
  cqm_path <- file.path(project_root, "data", "derived", "cell_qc_metrics.parquet")
  if (!file.exists(scores_path)) stop(sprintf("'%s' not found. Run `05_preprocessing_and_normalisation/03_calculate_program_scores.py` first.", scores_path))
  if (!file.exists(cqm_path)) stop(sprintf("'%s' not found. Run `04_quality_control/00_compute_cell_level_qc_metrics.py` first.", cqm_path))

  scores <- read_parquet(scores_path)
  cqm <- read_parquet(cqm_path)
  colnames(scores)[colnames(scores) == "__index_level_0__"] <- "cell_id"
  colnames(cqm)[colnames(cqm) == "__index_level_0__"] <- "cell_id"

  merged <- merge(
    scores,
    cqm[, c("cell_id", "section_id", "patient_id", "transcript_counts", "control_probe_ratio")],
    by = "cell_id"
  )

  set.seed(RNG_SEED)
  merged$row_id <- seq_len(nrow(merged))
  sampled_ids <- unlist(lapply(
    split(merged$row_id, merged$section_id),
    function(ids) sample(ids, min(SUBSAMPLE_PER_SECTION, length(ids)))
  ))
  df <- merged[sampled_ids, ]

  df$depth_z <- as.numeric(scale(log1p(df$transcript_counts)))
  df$control_burden_z <- as.numeric(scale(df$control_probe_ratio))
  df$patient_id <- factor(df$patient_id)
  df$section_id <- factor(df$section_id)

  cat(sprintf(
    "[INFO] Fitting mixed models on %d subsampled cells (%d sections, %d patients)...\n",
    nrow(df), length(unique(df$section_id)), length(unique(df$patient_id))
  ))

  results <- do.call(rbind, lapply(PROGRAMS, function(p) {
    cat(sprintf("[INFO] %s...\n", p))
    fit_and_decompose(df, p)
  }))

  output_dir <- file.path(project_root, "reports", "preprocess")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    results,
    file.path(project_root, "data", "derived", "variance_partition_summary.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  # Narrower, taller canvas: a word processor auto-scales a pasted
  # image to fit the page's text width (~6.5in), shrinking every
  # dimension including font size by that ratio. A narrower canvas
  # shrinks less, so the same in-image point size reads larger once
  # pasted.
  open_publication_pdf(file.path(output_dir, "variance_partition.pdf"), width = 9.8, height = 8.4)
  p <- ggplot(results, aes(x = program, y = fraction, fill = component)) +
    geom_col(position = "stack", width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = c(OKABE_ITO$blue, OKABE_ITO$orange, OKABE_ITO$bluish_green, OKABE_ITO$grey)) +
    labs(
      title = "Variance partition of preprocessing-stage\nprogramme scores",
      subtitle = "Diagnostic only: a technical QC check,\nnot a regression-adjustment step",
      x = NULL, y = "Fraction of variance", fill = "Component"
    ) +
    theme_publication()
  print(p)
  dev.off()

  n_unconverged <- sum(!results$converged[!duplicated(results$program)])
  cat("[OK] Variance partition summary:\n")
  print(results[, c("program", "component", "fraction")])
  cat(sprintf(
    "[OK] %d/%d model(s) reported a convergence warning. Wrote %s and data/derived/variance_partition_summary.tsv\n",
    n_unconverged, length(PROGRAMS), file.path(output_dir, "variance_partition.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
