#!/usr/bin/env Rscript
# `06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R` -- 05_resolve_myeloid_and_stromal_substates.R
#
# For each of 5 compartments (Myeloid, Dendritic_cell, Fibroblast,
# Endothelial, Perivascular_SmoothMuscle), gates cells to that compartment
# via `06_cell_type_annotation/02_score_major_lineages.py`'s argmax major-lineage call, then assigns the substate
# with the highest marker score among that compartment's own candidates
# (simple argmax within compartment -- unlike `06_cell_type_annotation/04_resolve_t_cell_substates.R`'s cross-cutting
# T-cell states, these substates are mutually exclusive by construction
# within their own compartment, so no priority-ordering scheme is needed
# here).
#
# Unlike `06_cell_type_annotation/04_resolve_t_cell_substates.R`, there is no external reference dataset for this
# compartment (GSE287301 is T cells only) -- no reference-transfer
# cross-check is possible or attempted here; this is panel-marker-evidence
# only, a documented scope difference, not an oversight.
#
# The actual marker-score computation runs in Python
# (_05_prepare_tme_substate_inputs.py, invoked below via system2()), not R,
# for the same reason as `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`, `06_cell_type_annotation/04_resolve_t_cell_substates.R` (no R HDF5/AnnData reader).
#
# Primary output: data/derived/tme_substates.parquet

suppressPackageStartupMessages({
  library(arrow)
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

# Compartment -> its own candidate substates, matching
# tme_substates.py::SUBSTATE_MARKERS exactly (column names are
# "{compartment}__{substate}_score").
COMPARTMENT_SUBSTATES <- list(
  Myeloid = c("Macrophage", "Monocyte"),
  Dendritic_cell = c("cDC", "pDC", "Mature_DC"),
  Fibroblast = c("Activated_CAF", "Resting_Fibroblast"),
  Endothelial = c("Blood_endothelial", "Lymphatic_endothelial"),
  Perivascular_SmoothMuscle = c("Pericyte", "Smooth_muscle")
)

# Substate marker sets are not all the same size (e.g. Mature_DC is a
# single gene, LAMP3, vs cDC's 5 genes) -- scanpy's score_genes does not
# average away noise the same way for a 1-gene list as it does for a
# multi-gene list, so raw scores across substates of different set sizes
# are not on a directly comparable scale. Confirmed on the data before
# fixing: Mature_DC's raw score had mean 0.774 vs cDC's 0.137 and pDC's
# (similarly small) mean, producing an implausible 75.7% "Mature_DC" call
# rate within the Dendritic_cell compartment (mature/migratory DCs are
# typically a tissue minority, not a majority, since they are actively
# migrating out to lymph nodes) -- the same class of scale-comparability
# problem, and the same z-score fix, as `06_cell_type_annotation/03_map_external_scrna_reference.py`'s cross-platform
# standardisation bug.
# Each substate's score is z-scored within its own compartment's cell
# population before comparing, not compared on its raw, un-standardised
# scale.
zscore <- function(x) {
  s <- sd(x)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x)) / s
}

assign_substate_within_compartment <- function(df, compartment, substates) {
  score_cols <- paste0(compartment, "__", substates, "_score")
  score_matrix <- as.matrix(df[, score_cols, drop = FALSE])
  score_matrix_z <- apply(score_matrix, 2, zscore)
  if (is.null(dim(score_matrix_z))) score_matrix_z <- matrix(score_matrix_z, ncol = length(substates))
  best_idx <- max.col(score_matrix_z, ties.method = "first")
  substates[best_idx]
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  helper_script <- file.path(project_root, "scripts", "06_cell_type_annotation", "_05_prepare_tme_substate_inputs.py")
  cat("[INFO] Running Python marker-score preparation helper...\n")
  exit_code <- system2("python3", c(shQuote(helper_script), "--project-root", shQuote(project_root)))
  if (exit_code != 0) {
    stop("Python marker-score preparation helper failed -- see its output above.")
  }

  inputs_path <- file.path(project_root, "data", "derived", "tme_substate_inputs.parquet")
  df <- read_parquet(inputs_path)
  colnames(df)[colnames(df) == "__index_level_0__"] <- "cell_id"

  df$tme_substate <- NA_character_
  for (compartment in names(COMPARTMENT_SUBSTATES)) {
    mask <- df$argmax_lineage == compartment
    if (sum(mask) == 0) next
    df$tme_substate[mask] <- assign_substate_within_compartment(
      df[mask, ], compartment, COMPARTMENT_SUBSTATES[[compartment]]
    )
  }

  gated <- df[!is.na(df$tme_substate), ]

  output_path <- file.path(project_root, "data", "derived", "tme_substates.parquet")
  write_parquet(gated[, c("cell_id", "argmax_lineage", "tme_substate")], output_path)

  cat("[INFO] Substate distribution per compartment:\n")
  for (compartment in names(COMPARTMENT_SUBSTATES)) {
    compartment_rows <- gated[gated$argmax_lineage == compartment, ]
    cat(sprintf("  %s (n=%d):\n", compartment, nrow(compartment_rows)))
    print(table(compartment_rows$tme_substate))
  }
  cat(sprintf(
    "[OK]   %d cell(s) total, %d gated across 5 compartments. Wrote %s\n",
    nrow(df), nrow(gated), output_path
  ))
}

if (sys.nframe() == 0) {
  main()
}
