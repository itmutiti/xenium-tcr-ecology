#!/usr/bin/env Rscript
# `06_cell_type_annotation/04_resolve_t_cell_substates.R` -- 04_resolve_t_cell_substates.R
#
# Assigns each T/NK-lineage cell (`06_cell_type_annotation/02_score_major_lineages.py` argmax) one discrete substate
# -- Treg, Cycling, Cytotoxic, Exhausted, CD4, CD8, or Ambiguous -- via a
# hierarchical priority rule over marker scores (`05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s cytotoxicity/
# exhaustion/proliferation scores, reused directly; a new Treg score; and
# direct CD4 vs CD8A expression). Cross-checks agreement against Phase
# 6.03's independently-derived reference-transfer labels wherever the two
# category systems overlap (Treg, Cycling, Exhausted, Cytotoxic).
#
# Priority order, most to least specific: Treg (FOXP3/IL2RA/CTLA4-positive
# AND CD4-dominant -- Tregs are conventionally CD4+, so CD4 dominance is
# required, not just a positive Treg score alone) > Cycling (a cell's
# expression profile is typically dominated by cell-cycle genes when
# actively proliferating, so this is checked before finer distinctions) >
# Exhausted > Cytotoxic > CD4/CD8 (by direct marker dominance) > Ambiguous
# (no marker score is clearly positive). This ordering is a documented
# judgment call -- not asserted to be the only defensible one -- since
# these programs are not mutually exclusive in T-cell biology (e.g. an
# exhausted cell can also be cycling).
#
# The actual marker-score computation runs in Python
# (_04_prepare_t_cell_substate_inputs.py, invoked below via system2()), not
# R, since this needs analysis_ready.h5ad and no R HDF5/AnnData reader is
# available in this project's environment (see `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`).
#
# Primary output: data/derived/t_cell_states.parquet

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

# Same "score_genes score > 0 means above the random background control
# gene set" convention already used implicitly across `05_preprocessing_and_normalisation/03_calculate_program_scores.py`, `06_cell_type_annotation/02_score_major_lineages.py`.
assign_t_cell_substate <- function(df) {
  is_treg <- df$treg_score > 0 & df$cd4_expr > df$cd8a_expr
  is_cycling <- df$proliferation_score > 0
  is_exhausted <- df$exhaustion_score > 0
  is_cytotoxic <- df$cytotoxicity_score > 0
  is_cd4 <- df$cd4_expr > df$cd8a_expr & df$cd4_expr > 0
  is_cd8 <- df$cd8a_expr > df$cd4_expr & df$cd8a_expr > 0

  state <- rep("Ambiguous", nrow(df))
  state[is_cd8] <- "CD8"
  state[is_cd4] <- "CD4"
  state[is_cytotoxic] <- "Cytotoxic"
  state[is_exhausted] <- "Exhausted"
  state[is_cycling] <- "Cycling"
  state[is_treg] <- "Treg"
  state
}

# Category systems overlap only partially between this script's states and
# `06_cell_type_annotation/03_map_external_scrna_reference.py`'s reference-transfer states -- CD4/CD8 (panel-only) and
# Naive_CM/Memory_TRM/MAIT_unconventional (reference-only) have no
# counterpart, so agreement is only meaningfully checked for the 4 states
# both systems name.
REFERENCE_STATE_MAP <- c(
  "Treg" = "Treg",
  "Cycling" = "Cycling",
  "Exhausted" = "Exhausted",
  "Cytotoxic" = "Cytotoxic_effector"
)

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  helper_script <- file.path(project_root, "scripts", "06_cell_type_annotation", "_04_prepare_t_cell_substate_inputs.py")
  cat("[INFO] Running Python marker-score preparation helper...\n")
  exit_code <- system2("python3", c(shQuote(helper_script), "--project-root", shQuote(project_root)))
  if (exit_code != 0) {
    stop("Python marker-score preparation helper failed -- see its output above.")
  }

  inputs_path <- file.path(project_root, "data", "derived", "t_cell_substate_inputs.parquet")
  df <- read_parquet(inputs_path)
  colnames(df)[colnames(df) == "__index_level_0__"] <- "cell_id"

  df$t_cell_state <- assign_t_cell_substate(df)

  gated <- df[df$argmax_lineage %in% c("T_cell", "NK_cell"), ]

  comparable <- gated[gated$t_cell_state %in% names(REFERENCE_STATE_MAP), ]
  expected_reference_state <- REFERENCE_STATE_MAP[comparable$t_cell_state]
  agreement <- mean(comparable$reference_predicted_state == expected_reference_state, na.rm = TRUE)

  # Written from `gated` (T/NK-lineage cells only), not the full `df`: the
  # hierarchical rule is computed for every cell (needed for the reference-
  # transfer agreement check above), but a T-cell substate call is only
  # meaningful for cells actually gated as T/NK lineage -- matching Phase
  # 6.05's tme_substates.parquet convention. An earlier version wrote the
  # full, ungated `df` here, which silently gave every cell in the dataset
  # (including e.g. Epithelial_Tumour cells) a "T-cell substate" label;
  # caught downstream in `06_cell_type_annotation/06_integrate_annotation_evidence.py` when the count of cells retaining a
  # substate call exceeded the total number of cells ever eligible for one.
  output_path <- file.path(project_root, "data", "derived", "t_cell_states.parquet")
  write_parquet(gated[, c("cell_id", "t_cell_state", "treg_score", "cd4_expr", "cd8a_expr",
                           "cytotoxicity_score", "exhaustion_score", "proliferation_score",
                           "argmax_lineage")], output_path)

  cat("[INFO] State distribution among gated T/NK-lineage cells:\n")
  print(table(gated$t_cell_state))
  msg <- paste0(
    "[OK]   %d cell(s) total, %d gated T/NK-lineage cell(s). Agreement with `06_cell_type_annotation/03_map_external_scrna_reference.py` reference ",
    "transfer (%d comparable-category cell(s)): %.3f. Wrote %s\n"
  )
  cat(sprintf(msg, nrow(df), nrow(gated), nrow(comparable), agreement, output_path))
}

if (sys.nframe() == 0) {
  main()
}
