#!/usr/bin/env Rscript
# `13_clone_ecology_confirmatory_models/02_fit_hierarchical_clone_models.R` -- 02_fit_hierarchical_clone_models.R
#
# Fits mixed/Bayesian models for any confirmed categorical structure,
# nested in patient and section -- explicitly conditional on Phase
# 11.05's prespecified structure-test result.
#
# Handling of the condition: `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s
# `q2_discrete_vs_continuous_structure_test` concluded continuous
# structure (frozen as taxonomy_version=v1_provisional, `13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py`) -- there is
# no confirmed categorical structure for this script to model. Fabricating
# a discrete clustering here to produce some output would directly
# contradict this project's prespecified confirmatory result. This
# script instead checks the frozen structure_type and, finding it
# continuous, writes an explicit "not applicable" record and report --
# the same "genuine gap, documented rather than hidden or faked"
# discipline already used for `06_cell_type_annotation/07_blinded_annotation_review.py`'s blinded
# panels.
#
# Primary output: reports/clone_ecology/hierarchical_models.pdf

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

#' Pure, testable: TRUE only if every analysis in the frozen
#' structure-test results confirms discrete/categorical structure.
check_categorical_structure_confirmed <- function(structure_test_results) {
  all(structure_test_results$final_decision == "discrete")
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  structure_test_path <- file.path(project_root, "data", "releases", "v1_clone_structure", "clone_structure_test_results.parquet")
  if (!file.exists(structure_test_path)) stop(sprintf("'%s' not found. Run `13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py` first.", structure_test_path))

  structure_test_results <- as.data.frame(read_parquet(structure_test_path))
  categorical_confirmed <- check_categorical_structure_confirmed(structure_test_results)

  output_dir_data <- file.path(project_root, "data", "derived")
  output_dir_reports <- file.path(project_root, "reports", "clone_ecology")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)

  if (!categorical_confirmed) {
    result <- data.frame(
      status = "not_applicable",
      reason = "`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s prespecified q2_discrete_vs_continuous_structure_test concluded continuous structure (both primary and sensitivity analyses); no categorical structure was confirmed for this script to model.",
      final_decisions = paste(structure_test_results$analysis, "=", structure_test_results$final_decision, collapse = "; ")
    )
    write_parquet(result, file.path(output_dir_data, "hierarchical_model_results.parquet"))

    # Narrower, taller canvas: a word processor auto-scales a pasted
    # image to fit the page's text width (~6.5in), shrinking every
    # dimension including font size by that ratio. A narrower canvas
    # shrinks less, so the same in-image point size reads larger once
    # pasted.
    open_publication_pdf(file.path(output_dir_reports, "hierarchical_models.pdf"), width = 9.4, height = 6.6)
    p <- ggplot() +
      xlim(0, 1) + ylim(0, 1) +
      annotate("text", x = 0.5, y = 0.72, label = "Hierarchical categorical-structure models", size = 8.3, fontface = "plain", family = "Liberation Sans") +
      annotate("text", x = 0.5, y = 0.54, label = "Not applicable", size = 9.1, fontface = "bold", colour = PUB_COLORS$sensitivity_analysis, family = "Liberation Sans") +
      annotate("text", x = 0.5, y = 0.34, label = "The prespecified discrete-vs-continuous structure test (Fig. 3) concluded\ncontinuous clone ecological structure; no categorical structure exists to model.", size = 6.2, family = "Liberation Sans") +
      annotate("text", x = 0.5, y = 0.14, label = "See Figure 3 (structure test) and Figure 4 (variance partition) for the confirmed structure.", size = 5.8, colour = "grey40", family = "Liberation Sans") +
      theme_void()
    print(p)
    dev.off()

    cat("[INFO] No confirmed categorical structure (`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`: continuous). This script is not applicable.\n")
    cat(sprintf("[OK]   Wrote %s, %s\n", file.path(output_dir_data, "hierarchical_model_results.parquet"), file.path(output_dir_reports, "hierarchical_models.pdf")))
    return(invisible(NULL))
  }

  stop("Categorical structure is confirmed but no discrete-model-fitting branch has been implemented yet -- this needs new code, not a placeholder, before it can run.")
}

if (sys.nframe() == 0) {
  main()
}
