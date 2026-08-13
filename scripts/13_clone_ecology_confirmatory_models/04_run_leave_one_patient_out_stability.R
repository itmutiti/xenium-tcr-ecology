#!/usr/bin/env Rscript
# `13_clone_ecology_confirmatory_models/04_run_leave_one_patient_out_stability.R` -- 04_run_leave_one_patient_out_stability.R
#
# Confirms `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`'s variance partition is not driven by any
# single patient -- refits the same nested model with each of the 10
# patients withheld in turn, checking whether all three variance
# components remain non-trivial in every fold.
#
# Scope note: "any categorical structure" (this script's scaffold
# wording) does not apply here -- `13_clone_ecology_confirmatory_models/02_fit_hierarchical_clone_models.R` already found, and
# documented, that `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R` confirmed continuous, not
# categorical, structure; there is no categorical-structure LOPO check to run
# alongside the variance-partition one.
#
# Primary output: reports/clone_ecology/lopo_stability.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(lme4)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

# A defensible, round threshold for "non-trivial": a component
# accounting for less than 5% of total variance is arguably negligible;
# the full-data point estimates (patient 29.3%, identity 20.2%, context
# 50.4%, `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`) are all far above this, so the question this
# threshold operationalises is whether removing any one patient can push
# a component down near/below it.
NON_TRIVIAL_THRESHOLD <- 0.05

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

#' Pure, testable: TRUE if every one of the three proportions clears
#' `threshold`.
all_components_non_trivial <- function(proportions, threshold = NON_TRIVIAL_THRESHOLD) {
  all(c(proportions$patient_proportion, proportions$identity_proportion, proportions$context_proportion) >= threshold)
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  structure_path <- file.path(project_root, "data", "releases", "v1_clone_structure", "clone_ecological_structure.parquet")
  if (!file.exists(structure_path)) stop(sprintf("'%s' not found. Run `13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py` first.", structure_path))

  data <- as.data.frame(read_parquet(structure_path))
  patients <- sort(unique(data$patient_id))
  cat(sprintf("[INFO] Data: n=%d clone-section rows, %d patients.\n", nrow(data), length(patients)))

  rows <- list()
  for (p in patients) {
    held_out_data <- data[data$patient_id != p, ]
    model <- tryCatch(
      lme4::lmer(ecological_structure_score ~ 1 + (1 | patient_id/clone_id), data = held_out_data, REML = TRUE),
      error = function(e) NULL, warning = function(w) NULL
    )
    if (is.null(model)) {
      rows[[p]] <- data.frame(withheld_patient = p, patient_proportion = NA, identity_proportion = NA, context_proportion = NA, all_non_trivial = NA, converged = FALSE)
      next
    }
    vc <- extract_variance_partition(model)
    prop <- compute_variance_proportions(vc$patient_var, vc$identity_var, vc$context_var)
    rows[[p]] <- data.frame(
      withheld_patient = p,
      patient_proportion = prop$patient_proportion,
      identity_proportion = prop$identity_proportion,
      context_proportion = prop$context_proportion,
      all_non_trivial = all_components_non_trivial(prop),
      converged = TRUE
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  n_stable_folds <- sum(result$all_non_trivial, na.rm = TRUE)
  n_folds <- nrow(result)

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "variance_partition_lopo_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "clone_ecology")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Narrower, taller canvas: a word processor auto-scales a pasted
  # image to fit the page's text width (~6.5in), shrinking every
  # dimension including font size by that ratio. A narrower canvas
  # shrinks less, so the same in-image point size reads larger once
  # pasted.
  open_publication_pdf(file.path(output_dir_reports, "lopo_stability.pdf"), width = 9.6, height = 7.6)

  plot_data <- rbind(
    data.frame(withheld_patient = result$withheld_patient, proportion = result$patient_proportion, component = "patient"),
    data.frame(withheld_patient = result$withheld_patient, proportion = result$identity_proportion, component = "identity"),
    data.frame(withheld_patient = result$withheld_patient, proportion = result$context_proportion, component = "context")
  )
  plot_data$component <- factor(plot_data$component, levels = c("patient", "identity", "context"), labels = c("Patient", "Clonal identity", "Spatial context"))
  p1 <- ggplot(plot_data, aes(x = withheld_patient, y = proportion, fill = component)) +
    geom_col(position = "stack", width = 0.7) +
    geom_hline(yintercept = NON_TRIVIAL_THRESHOLD, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    scale_fill_manual(values = c("Patient" = OKABE_ITO$blue, "Clonal identity" = OKABE_ITO$orange, "Spatial context" = OKABE_ITO$bluish_green)) +
    labs(
      subtitle = sprintf("%d/%d folds retain all three components >= %.0f%% of variance", n_stable_folds, n_folds, NON_TRIVIAL_THRESHOLD * 100),
      x = "Withheld patient", y = "Proportion of total variance", fill = NULL
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1)

  dev.off()

  cat(sprintf("[INFO] %d/%d LOPO folds retain all three components >= %.0f%% non-trivial threshold.\n", n_stable_folds, n_folds, NON_TRIVIAL_THRESHOLD * 100))
  for (i in seq_len(nrow(result))) {
    cat(sprintf(
      "[INFO]   withheld=%-6s patient=%.3f identity=%.3f context=%.3f non_trivial=%s\n",
      result$withheld_patient[i], result$patient_proportion[i], result$identity_proportion[i], result$context_proportion[i], result$all_non_trivial[i]
    ))
  }
  cat(sprintf(
    "[OK]   LOPO stability check complete. Wrote %s, %s\n",
    file.path(output_dir_data, "variance_partition_lopo_results.parquet"),
    file.path(output_dir_reports, "lopo_stability.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
