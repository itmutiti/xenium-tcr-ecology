#!/usr/bin/env Rscript
# `04_quality_control/08_assess_replicate_concordance.R` -- 08_assess_replicate_concordance.R
#
# Quantifies technical replicate agreement between the two Xenium runs for
# each of the 7 replicated patients (P09, P12, P13, P17, P19, P20, P28,
# confirmed against the data -- matches McCord et al. 2026's count of
# "seven patient samples ... profiled in duplicate"), across four axes:
# cell/count-level agreement, pseudobulk gene expression correlation,
# patient-specific TCR/CDR3 probe detection concordance (the paper's own
# metric, "1 of 7 pairs discordant" -- McCord et al. 2026, Sci Immunol
# 11:eaec3133), and a self-contained spatial autocorrelation (Moran's I)
# comparison.
#
# Deliberately not a "day nested within patient" mixed-effects model: with
# exactly 2 sections (1 day-pair) per patient and no further nesting to
# partition, there is no variance-component structure a mixed model would
# add over a direct paired comparison -- that random-effects framing
# belongs to later phases that pool replicate and non-replicate sections
# together (e.g. Clone Ecology Confirmatory Models's variance-partition models), not this per-pair
# concordance report.
#
# Primary output: reports/qc/replicate_concordance.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(Matrix)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

# Mirrors src/xenium_tcr_ecology/infra/paths.py::find_project_root (see
# `04_quality_control/06_define_qc_thresholds_hierarchically.R` for the first R port of this convention).
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

# Python's pandas writes boolean columns as the literal strings "True"/
# "False" -- read.delim() on a TSV does not recognise these as logical (only
# "TRUE"/"FALSE" are auto-detected), silently leaving a character column;
# using that directly as a `df[mask, ]` row index is interpreted by R as
# name-matching, not as a logical mask, and silently returns all-NA rows.
# arrow::read_parquet() has its own version of the same hazard: it returns
# these columns as a *factor* (confirmed against the data), not logical.
# as_bool() normalises either representation. Confirmed as the root cause of
# an initial NA result in this script's own Moran's I computation during
# development, section 9.
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  as.character(x) %in% c("True", "TRUE", "T", "true")
}

MORAN_K_NEIGHBOURS <- 6
SPATIAL_SUBSAMPLE_SIZE <- 2000
# Same robust-outlier convention as `04_quality_control/02_detect_spatial_qc_artifacts.py`, `04_quality_control/06_define_qc_thresholds_hierarchically.R` (Iglewicz & Hoaglin, 1993).
MAD_FLAG_THRESHOLD <- 3.5
# Confirmed against the panel's gene names.
# The date prefix is 6 digits optionally followed by a single batch letter
# (e.g. "231004B") -- an earlier version of this pattern omitted the
# optional letter and silently missed 17 CDR3 probes (the entire
# "231004B" batch, present in P12 and P28's panels), undercounting those two
# patients' probe sets in this script's first run. Fixed and re-run, section 12.
CDR3_PROBE_PATTERN <- "^[0-9]{6}[A-Z]?_[A-Z]+_TR[AB]$"

# --- Pure, testable helper functions -----------------------------------

moran_i <- function(coords, values, k = MORAN_K_NEIGHBOURS) {
  n <- nrow(coords)
  d <- as.matrix(dist(coords))
  w <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nn <- order(d[i, ])[2:(k + 1)]
    w[i, nn] <- 1
  }
  w <- w / rowSums(w)
  dev <- values - mean(values)
  sum(w * outer(dev, dev)) / sum(dev^2)
}

jaccard_index <- function(set_a, set_b) {
  u <- union(set_a, set_b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(set_a, set_b)) / length(u)
}

modified_z_scores <- function(x) {
  med <- median(x)
  mad_val <- median(abs(x - med))
  if (mad_val == 0) return(rep(0, length(x)))
  0.6745 * (x - med) / mad_val
}

hull_area <- function(x, y) {
  hpts <- chull(x, y)
  hpts <- c(hpts, hpts[1])
  xh <- x[hpts]
  yh <- y[hpts]
  0.5 * abs(sum(xh[-length(xh)] * yh[-1] - xh[-1] * yh[-length(yh)]))
}

# --- Data access ---------------------------------------------------------

read_exclusion_log <- function(project_root) {
  path <- file.path(project_root, "data", "derived", "exclusion_log.tsv")
  if (!file.exists(path)) stop(sprintf("'%s' not found. Run `04_quality_control/07_apply_qc_filters_with_audit_trail.py` first.", path))
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE)
  df$qc_pass <- as_bool(df$qc_pass)
  df
}

read_qc_pass_obs <- function(project_root, section_id, exclusion_log) {
  dir <- file.path(project_root, "data", "objects", "r_exports", section_id)
  obs <- read_parquet(file.path(dir, "obs_metadata.parquet"))
  excl_section <- exclusion_log[exclusion_log$section_id == section_id, ]
  local_id <- sub(paste0("^", section_id, "_"), "", excl_section$cell_id)
  m <- match(obs$cell_id, local_id)
  obs$qc_pass <- excl_section$qc_pass[m]
  obs[!is.na(obs$qc_pass) & obs$qc_pass, ]
}

read_section_matrix <- function(project_root, section_id) {
  dir <- file.path(project_root, "data", "objects", "r_exports", section_id)
  mat <- readMM(file.path(dir, "matrix.mtx.gz"))
  barcodes <- readLines(gzfile(file.path(dir, "barcodes.tsv.gz")))
  genes <- read_parquet(file.path(dir, "var_metadata.parquet"))$gene_ids
  rownames(mat) <- genes
  colnames(mat) <- barcodes
  mat
}

pseudobulk_for_cells <- function(mat, cell_ids) {
  keep <- intersect(colnames(mat), cell_ids)
  pb <- Matrix::rowSums(mat[, keep, drop = FALSE])
  pb
}

discover_replicate_pairs <- function(project_root) {
  r_exports_root <- file.path(project_root, "data", "objects", "r_exports")
  sections <- list.dirs(r_exports_root, full.names = FALSE, recursive = FALSE)
  meta <- do.call(rbind, lapply(sections, function(s) {
    df <- read_parquet(
      file.path(r_exports_root, s, "obs_metadata.parquet"),
      col_select = c("patient_id", "run_number", "is_technical_replicate")
    )
    data.frame(
      section_id = s,
      patient_id = as.character(df$patient_id[1]),
      run_number = as.integer(as.character(df$run_number[1])),
      is_technical_replicate = as_bool(df$is_technical_replicate[1]),
      stringsAsFactors = FALSE
    )
  }))
  replicate_meta <- meta[meta$is_technical_replicate, ]
  patients <- unique(replicate_meta$patient_id)
  pairs <- do.call(rbind, lapply(patients, function(p) {
    rows <- replicate_meta[replicate_meta$patient_id == p, ]
    rows <- rows[order(rows$run_number), ]
    if (nrow(rows) != 2) {
      stop(sprintf(
        "Patient '%s' is marked as a technical replicate but has %d section(s), not 2 -- cannot form a pair.",
        p, nrow(rows)
      ))
    }
    data.frame(patient_id = p, section1 = rows$section_id[1], section2 = rows$section_id[2], stringsAsFactors = FALSE)
  }))
  pairs
}

# --- Core per-pair assessment --------------------------------------------

assess_pair <- function(project_root, patient_id, section1, section2, exclusion_log) {
  obs1 <- read_qc_pass_obs(project_root, section1, exclusion_log)
  obs2 <- read_qc_pass_obs(project_root, section2, exclusion_log)

  mat1 <- read_section_matrix(project_root, section1)
  mat2 <- read_section_matrix(project_root, section2)

  pb1 <- pseudobulk_for_cells(mat1, obs1$cell_id)
  pb2 <- pseudobulk_for_cells(mat2, obs2$cell_id)

  shared_genes <- intersect(names(pb1), names(pb2))
  pseudobulk_r <- cor(log1p(pb1[shared_genes]), log1p(pb2[shared_genes]), method = "pearson")

  cdr3_genes <- shared_genes[grepl(CDR3_PROBE_PATTERN, shared_genes)]
  detected1 <- cdr3_genes[pb1[cdr3_genes] > 0]
  detected2 <- cdr3_genes[pb2[cdr3_genes] > 0]
  cdr3_jaccard <- jaccard_index(detected1, detected2)

  set.seed(RNG_SEED)
  idx1 <- sample(seq_len(nrow(obs1)), min(SPATIAL_SUBSAMPLE_SIZE, nrow(obs1)))
  idx2 <- sample(seq_len(nrow(obs2)), min(SPATIAL_SUBSAMPLE_SIZE, nrow(obs2)))
  moran1 <- moran_i(as.matrix(obs1[idx1, c("x_centroid", "y_centroid")]), obs1$transcript_counts[idx1])
  moran2 <- moran_i(as.matrix(obs2[idx2, c("x_centroid", "y_centroid")]), obs2$transcript_counts[idx2])

  area1 <- hull_area(obs1$x_centroid, obs1$y_centroid)
  area2 <- hull_area(obs2$x_centroid, obs2$y_centroid)

  data.frame(
    patient_id = patient_id,
    section1 = section1,
    section2 = section2,
    n_cells_1 = nrow(obs1),
    n_cells_2 = nrow(obs2),
    cell_density_1_per_mm2 = nrow(obs1) / (area1 / 1e6),
    cell_density_2_per_mm2 = nrow(obs2) / (area2 / 1e6),
    median_counts_1 = median(obs1$transcript_counts),
    median_counts_2 = median(obs2$transcript_counts),
    pseudobulk_pearson_r = pseudobulk_r,
    n_cdr3_probes_shared = length(cdr3_genes),
    n_cdr3_detected_1 = length(detected1),
    n_cdr3_detected_2 = length(detected2),
    cdr3_probe_jaccard = cdr3_jaccard,
    moran_i_1 = moran1,
    moran_i_2 = moran2,
    moran_i_abs_diff = abs(moran1 - moran2),
    stringsAsFactors = FALSE
  )
}

# --- Orchestration --------------------------------------------------------

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))
  exclusion_log <- read_exclusion_log(project_root)
  pairs <- discover_replicate_pairs(project_root)

  cat(sprintf("[INFO] Found %d replicate pair(s): %s\n", nrow(pairs), paste(pairs$patient_id, collapse = ", ")))

  results <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    row <- pairs[i, ]
    cat(sprintf("[INFO] Assessing pair %s (%s vs %s)...\n", row$patient_id, row$section1, row$section2))
    assess_pair(project_root, row$patient_id, row$section1, row$section2, exclusion_log)
  }))

  # Primary concordance metric for outlier flagging: pseudobulk expression
  # correlation -- the most standard, directly comparable replicate-QC
  # statistic across all 7 pairs (unlike Moran's I or density, which are
  # tissue-content-dependent even between genuine biological replicates).
  results$pseudobulk_r_mad_z <- modified_z_scores(results$pseudobulk_pearson_r)
  results$flagged_discordant <- results$pseudobulk_r_mad_z < -MAD_FLAG_THRESHOLD

  output_dir <- file.path(project_root, "reports", "qc")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    results,
    file.path(project_root, "data", "derived", "replicate_concordance.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  # Mixed 2x2 grid (A, B on top; C on bottom-left; bottom-right left
  # blank by compose_panels()'s row-major fill for n=3), not a full
  # vertical stack or full horizontal row: a word processor auto-scales
  # a pasted image to fit the page's text width (~6.5in), shrinking
  # every dimension including font size by that ratio, so an
  # all-horizontal 3-in-a-row image shrinks font size the most; an
  # all-vertical 3-stack becomes excessively long. The 2x2 compromise
  # keeps the width (and therefore the shrink ratio) moderate while
  # keeping each panel roughly as tall as it needs to be.
  open_publication_pdf(file.path(output_dir, "replicate_concordance.pdf"), width = 16.0, height = 13.5)

  p1 <- ggplot(results, aes(x = reorder(patient_id, pseudobulk_pearson_r), y = pseudobulk_pearson_r, fill = flagged_discordant)) +
    geom_col(width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = PUB_COLORS$not_significant, "TRUE" = PUB_COLORS$sensitivity_analysis), guide = "none") +
    labs(
      subtitle = sprintf("Pseudobulk concordance (vermillion = MAD outlier, |z|>%.1f)", MAD_FLAG_THRESHOLD),
      x = "Patient", y = "Pearson r (run1 vs run2)"
    ) +
    theme_publication()

  p2 <- ggplot(results, aes(x = reorder(patient_id, -cdr3_probe_jaccard), y = cdr3_probe_jaccard)) +
    geom_col(fill = PUB_COLORS$primary_analysis, width = 0.65) +
    coord_flip() +
    labs(
      subtitle = "CDR3-probe detection concordance (Jaccard index)",
      x = NULL, y = "Jaccard index"
    ) +
    theme_publication()

  p3 <- ggplot(results, aes(x = moran_i_1, y = moran_i_2, label = patient_id)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    geom_point(size = 3.6, colour = PUB_COLORS$primary_analysis) +
    geom_text(vjust = -1, size = 5.7, family = "Liberation Sans") +
    # Extra headroom above the highest point: its text label (vjust=-1,
    # above the point) was otherwise clipped by the top of the panel.
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(
      subtitle = sprintf("Spatial autocorrelation concordance (k = %d, %d-cell subsample)", MORAN_K_NEIGHBOURS, SPATIAL_SUBSAMPLE_SIZE),
      x = "Moran's I (run1)", y = "Moran's I (run2)"
    ) +
    theme_publication()

  compose_panels(list(p1, p2, p3), ncol = 2)

  dev.off()

  cat(sprintf(
    "[OK]   %d pair(s) assessed. %d flagged as discordant on pseudobulk expression: %s. Wrote %s and data/derived/replicate_concordance.tsv\n",
    nrow(results), sum(results$flagged_discordant),
    paste(results$patient_id[results$flagged_discordant], collapse = ", "),
    file.path(output_dir, "replicate_concordance.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
