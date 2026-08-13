#!/usr/bin/env Rscript
# `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R` -- 05_test_discrete_vs_continuous_structure.R
#
# The prespecified `q2_discrete_vs_continuous_structure_test`
# confirmatory analysis (governance/analysis_registry.tsv):
# "Clone ecological structure is better described by a continuous latent
# axis than by discrete categorical clusters" is prespecified as the
# null to beat, avoiding restating the source paper's illustrative,
# non-quantitative "spatial fate" claims on a handful of hand-picked
# clones (McCord et al. 2026, Sci Immunol 11:eaec3133). Discrete
# taxonomy is reported only
# if it clears this prespecified burden of proof; a tie or inconclusive
# result falls back to the continuous null.
#
# No specialised clustering/factor-analysis R package is installed in
# this environment, so every method below is implemented from first-
# principles statistics using only base R's `stats` package (`kmeans`,
# `factanal`, `prcomp`, `solve`, `determinant`).
#
# Primary endpoint (per the registry): a model-comparison statistic
# (BIC) and a gap statistic, evaluated independently -- discrete
# structure is adopted only if both independent methods agree.
#
# Scope note: this script produces the test result and report only.
# Freezing `taxonomy_version` in `governance/taxonomy_version_log.tsv` is
# explicitly `11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py`'s primary output, not this script's -- kept
# separate per the division of labour between "test structure" (11.05)
# and "freeze taxonomy version" (11.07).
#
# Primary output: reports/clone_ecology/structure_test.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

K_MAX <- 6
N_GAP_REFERENCE <- 50
DELTA_BIC_THRESHOLD <- 10 # Kass & Raftery (1995): "very strong evidence"

# Full-coverage (100% of 261 clone-section rows) feature set for the
# primary confirmatory analysis -- deliberately excludes the rarefaction-
# gated spatial-pattern descriptors (`11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`: 56-58% coverage) and
# barrier metrics (`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`: 58% coverage) from the primary test to
# maximise statistical power on the full population; those are instead
# included in a sensitivity analysis on the smaller complete-case
# subset. One state fraction (`ambiguous_fraction`) and one raw adjacency
# (`malignant_adjacency`, redundant with the opportunity-normalised
# `engagement_ratio`) are deliberately dropped to avoid a singular
# covariance matrix from perfectly collinear composition fractions.
PRIMARY_FEATURES <- c(
  "cytotoxic_fraction", "exhausted_fraction", "cycling_fraction", "treg_fraction", "cd4_fraction", "cd8_fraction",
  "shannon_entropy_bits", "engagement_ratio", "dc_engagement_ratio", "macrophage_engagement_ratio",
  "antigen_presentation_score_excess"
)
SENSITIVITY_EXTRA_FEATURES <- c(
  "dispersion_rarefied", "clark_evans_index_rarefied", "domain_richness_rarefied",
  "fibroblast_barrier_fraction", "vascular_barrier_fraction", "suppressive_myeloid_barrier_fraction"
)

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

#' Multivariate normal log-likelihood of mean-centred data under
#' N(0, sigma) -- base-R only (no mvtnorm dependency).
compute_mvn_loglik <- function(data_centered, sigma) {
  p <- ncol(data_centered)
  sigma_inv <- solve(sigma)
  log_det <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
  quad_form <- rowSums((data_centered %*% sigma_inv) * data_centered)
  sum(-0.5 * (p * log(2 * pi) + log_det + quad_form))
}

#' BIC of a k-means-approximated spherical Gaussian mixture at K=k, with
#' per-cluster variance (k=1 is a single Gaussian, i.e. "no discrete
#' substructure"). Pure, testable given a seeded k-means fit.
#'
#' A shared (pooled) variance was tried first and rejected: checked
#' directly against synthetic pure-noise data (n=100, p=2) before relying
#' on it, and it kept preferring higher K even for a single unstructured
#' Gaussian cloud (BIC fell monotonically K=1..5) -- because sigma^2
#' shrinks with K regardless of whether the extra clusters reflect
#' structure, and one shared-variance parameter is too cheap a penalty to
#' offset that. Per-cluster variance (each cluster's within-cluster
#' SS / (n_k * p)) is the more standard formulation (matching e.g.
#' Pelleg & Moore 2000's X-means BIC) and a larger parameter penalty
#' (k*p + k + (k-1) instead of k*p + 1 + (k-1)) that behaves far more
#' conservatively on the same synthetic checks.
compute_discrete_bic <- function(data, k, seed = RNG_SEED) {
  n <- nrow(data)
  p <- ncol(data)
  if (k == 1) {
    centre <- colMeans(data)
    ss <- sum(scale(data, center = centre, scale = FALSE)^2)
    sigma2 <- ss / (n * p)
    log_lik <- -0.5 * n * p * log(2 * pi * sigma2) - n * p / 2
    n_params <- p + 1 # one mean vector + one variance
  } else {
    set.seed(seed)
    km <- kmeans(data, centers = k, nstart = 25, iter.max = 100)
    log_lik <- 0
    for (cluster_id in seq_len(k)) {
      n_k <- sum(km$cluster == cluster_id)
      if (n_k < 2) next # too few points in this cluster to estimate a variance; contributes 0, not fabricated
      sigma2_k <- km$withinss[cluster_id] / (n_k * p)
      if (sigma2_k <= 0) next
      log_lik <- log_lik - 0.5 * n_k * p * log(2 * pi * sigma2_k) - n_k * p / 2
    }
    n_params <- k * p + k + (k - 1) # means + per-cluster variances + mixing proportions
  }
  -2 * log_lik + n_params * log(n)
}

#' BIC of a 1-factor continuous latent-axis model (`factanal`, base R).
#' Pure, testable.
compute_continuous_bic <- function(data, n_factors = 1) {
  n <- nrow(data)
  p <- ncol(data)
  data_scaled <- scale(data)
  fa <- factanal(data_scaled, factors = n_factors)
  loadings_matrix <- matrix(fa$loadings, nrow = p, ncol = n_factors)
  sigma_hat <- loadings_matrix %*% t(loadings_matrix) + diag(fa$uniquenesses)
  log_lik <- compute_mvn_loglik(data_scaled, sigma_hat)
  n_params <- p * (p + 1) / 2 - fa$dof
  -2 * log_lik + n_params * log(n)
}

#' Tibshirani et al. (2001) gap statistic, PCA-aligned reference
#' sampling. Pure, testable given a seeded reference draw.
compute_gap_statistic <- function(data, k_max = K_MAX, b_reference = N_GAP_REFERENCE, seed = RNG_SEED) {
  n <- nrow(data)
  p <- ncol(data)

  compute_log_wk <- function(x, k) {
    if (k == 1) {
      centre <- colMeans(x)
      return(log(sum(scale(x, center = centre, scale = FALSE)^2)))
    }
    km <- kmeans(x, centers = k, nstart = 10, iter.max = 100)
    log(km$tot.withinss)
  }

  observed_log_wk <- sapply(seq_len(k_max), function(k) compute_log_wk(data, k))

  pca <- prcomp(data, center = TRUE, scale. = FALSE)
  mins <- apply(pca$x, 2, min)
  maxs <- apply(pca$x, 2, max)

  ref_log_wk <- matrix(NA_real_, nrow = b_reference, ncol = k_max)
  for (b in seq_len(b_reference)) {
    set.seed(seed + b)
    ref_rotated <- sapply(seq_len(p), function(j) stats::runif(n, min = mins[j], max = maxs[j]))
    ref_data <- sweep(ref_rotated %*% t(pca$rotation), 2, pca$center, "+")
    for (k in seq_len(k_max)) {
      ref_log_wk[b, k] <- compute_log_wk(ref_data, k)
    }
  }

  gap <- colMeans(ref_log_wk) - observed_log_wk
  sd_k <- apply(ref_log_wk, 2, sd) * sqrt(1 + 1 / b_reference)
  data.frame(k = seq_len(k_max), gap = gap, sd = sd_k, observed_log_wk = observed_log_wk)
}

#' Tibshirani's standard selection rule: smallest k with
#' Gap(k) >= Gap(k+1) - sd(k+1). Pure, testable.
select_k_by_gap <- function(gap_df) {
  k_max <- nrow(gap_df)
  for (k in seq_len(k_max - 1)) {
    if (gap_df$gap[k] >= gap_df$gap[k + 1] - gap_df$sd[k + 1]) {
      return(k)
    }
  }
  k_max
}

#' Decision rule: discrete structure is adopted only if both the BIC
#' comparison and the gap statistic independently favour it -- the
#' prespecified continuous null is retained on any disagreement or
#' inconclusive result.
classify_structure <- function(delta_bic, gap_selected_k, delta_bic_threshold = DELTA_BIC_THRESHOLD) {
  bic_favors_discrete <- delta_bic > delta_bic_threshold
  gap_favors_discrete <- gap_selected_k >= 2
  if (bic_favors_discrete && gap_favors_discrete) {
    return("discrete")
  }
  "continuous"
}

run_structure_test <- function(data, features, label) {
  feature_matrix <- as.matrix(data[, features])
  n <- nrow(feature_matrix)

  discrete_bic <- sapply(1:K_MAX, function(k) compute_discrete_bic(feature_matrix, k))
  best_k <- which.min(discrete_bic)
  best_discrete_bic <- discrete_bic[best_k]
  continuous_bic <- compute_continuous_bic(feature_matrix, n_factors = 1)
  delta_bic <- continuous_bic - best_discrete_bic

  gap_df <- compute_gap_statistic(feature_matrix, k_max = K_MAX)
  gap_selected_k <- select_k_by_gap(gap_df)

  decision <- classify_structure(delta_bic, gap_selected_k)

  list(
    label = label,
    n = n,
    n_features = length(features),
    discrete_bic_by_k = discrete_bic,
    best_k = best_k,
    best_discrete_bic = best_discrete_bic,
    continuous_bic = continuous_bic,
    delta_bic = delta_bic,
    gap_df = gap_df,
    gap_selected_k = gap_selected_k,
    decision = decision
  )
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())

  paths <- list(
    spatial = file.path(project_root, "data", "derived", "clone_spatial_descriptors.parquet"),
    state = file.path(project_root, "data", "derived", "clone_state_composition.parquet"),
    engagement = file.path(project_root, "data", "derived", "clone_tumour_engagement.parquet"),
    apc = file.path(project_root, "data", "derived", "clone_apc_support.parquet"),
    barrier = file.path(project_root, "data", "derived", "clone_barrier_metrics.parquet")
  )
  for (p in paths) {
    if (!file.exists(p)) stop(sprintf("'%s' not found. Run `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`-11.04 first.", p))
  }

  spatial <- as.data.frame(read_parquet(paths$spatial))
  state <- as.data.frame(read_parquet(paths$state))
  engagement <- as.data.frame(read_parquet(paths$engagement))
  apc <- as.data.frame(read_parquet(paths$apc))
  barrier <- as.data.frame(read_parquet(paths$barrier))

  key <- c("clone_id", "section_id", "patient_id", "n_cells")
  merged <- Reduce(function(x, y) merge(x, y, by = key), list(spatial, state, engagement, apc, barrier))
  if (nrow(merged) != nrow(spatial)) {
    stop(sprintf("Merge changed row count (%d -> %d) -- `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`-11.04 populations are expected to match exactly.", nrow(spatial), nrow(merged)))
  }

  primary_data <- merged[stats::complete.cases(merged[, PRIMARY_FEATURES]), ]
  sensitivity_features <- c(PRIMARY_FEATURES, SENSITIVITY_EXTRA_FEATURES)
  sensitivity_data <- merged[stats::complete.cases(merged[, sensitivity_features]), ]

  cat(sprintf("[INFO] Primary analysis: n=%d clone-sections, %d features (full coverage).\n", nrow(primary_data), length(PRIMARY_FEATURES)))
  cat(sprintf("[INFO] Sensitivity analysis: n=%d clone-sections, %d features (complete-case, includes rarefied + barrier descriptors).\n", nrow(sensitivity_data), length(sensitivity_features)))

  primary_result <- run_structure_test(primary_data, PRIMARY_FEATURES, "primary")
  sensitivity_result <- run_structure_test(sensitivity_data, sensitivity_features, "sensitivity")

  for (r in list(primary_result, sensitivity_result)) {
    cat(sprintf(
      "[INFO] [%s] n=%d, best discrete K=%d (BIC=%.1f), continuous BIC=%.1f, delta_BIC=%.1f, gap-selected K=%d -> decision: %s\n",
      r$label, r$n, r$best_k, r$best_discrete_bic, r$continuous_bic, r$delta_bic, r$gap_selected_k, r$decision
    ))
  }

  final_decision <- if (primary_result$decision == "discrete" && sensitivity_result$decision == "discrete") "discrete" else "continuous"
  cat(sprintf("[INFO] Final decision (requires both primary and sensitivity to agree on 'discrete'): %s\n", final_decision))

  results_table <- data.frame(
    analysis = c("primary", "sensitivity"),
    n = c(primary_result$n, sensitivity_result$n),
    n_features = c(primary_result$n_features, sensitivity_result$n_features),
    best_discrete_k = c(primary_result$best_k, sensitivity_result$best_k),
    best_discrete_bic = c(primary_result$best_discrete_bic, sensitivity_result$best_discrete_bic),
    continuous_bic = c(primary_result$continuous_bic, sensitivity_result$continuous_bic),
    delta_bic = c(primary_result$delta_bic, sensitivity_result$delta_bic),
    gap_selected_k = c(primary_result$gap_selected_k, sensitivity_result$gap_selected_k),
    decision = c(primary_result$decision, sensitivity_result$decision),
    final_decision = final_decision
  )

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(results_table, file.path(output_dir_data, "clone_structure_test_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "clone_ecology")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)

  # Shared publication theme/colour convention, matching the Python-side
  # `xenium_tcr_ecology.viz.style` module used by every other main
  # figure (Okabe-Ito colourblind-safe palette, open axes, no gridlines,
  # Liberation Sans/Arial typography, neutral titles only).
  COL_CURVE <- "#0072B2" # Okabe-Ito blue: the criterion curve (BIC or gap)
  COL_SELECTED <- "#D55E00" # Okabe-Ito vermillion: the criterion's chosen K
  COL_REFERENCE <- "#999999" # Okabe-Ito grey: reference value

  # Typography/line-weight scale raised twice 2026-07-19 (base_size
  # 9.5->14->19, matching the equivalent increases in the shared
  # style.py/theme.R scale): first after feedback that the original
  # scale read as undersized in print, then again after feedback that
  # text was still too small once pasted into a word processor and
  # auto-scaled to fit a page width -- legibility after that scale-down
  # depends on the font-to-canvas-width ratio, so font size was raised
  # well ahead of the (unchanged) canvas size below.
  theme_publication <- function() {
    theme_minimal(base_family = "Liberation Sans", base_size = 22) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(colour = "grey30", linewidth = 0.6),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.6),
        axis.ticks.length = unit(4, "pt"),
        plot.title = element_text(size = 22, face = "plain", hjust = 0),
        plot.margin = margin(12, 14, 8, 8)
      )
  }

  # Feature-set context is now embedded directly in each panel's own
  # title (subtitle line) rather than via a shared rotated row label --
  # the vertical single-column layout below has no row grouping left
  # for a row label to annotate.
  panel_plot <- function(r, kind, feature_set_label) {
    if (kind == "bic") {
      bic_df <- data.frame(k = 1:K_MAX, bic = r$discrete_bic_by_k)
      p <- ggplot(bic_df, aes(x = k, y = bic)) +
        geom_hline(yintercept = r$continuous_bic, linetype = "dashed", colour = COL_REFERENCE, linewidth = 0.7) +
        geom_line(colour = COL_CURVE, linewidth = 1.0) +
        geom_point(colour = COL_CURVE, size = 2.6) +
        geom_point(data = bic_df[r$best_k, ], aes(x = k, y = bic), colour = COL_SELECTED, size = 4.2) +
        scale_x_continuous(breaks = 1:K_MAX) +
        labs(x = "Discrete clusters (K)", y = "BIC", title = sprintf("%s: BIC", feature_set_label)) +
        theme_publication()
    } else {
      # Error bars (+/-1 SD) are secondary context around the gap-statistic
      # curve, not the primary signal -- lightened so the curve and its
      # selected-K point remain the dominant visual element.
      p <- ggplot(r$gap_df, aes(x = k, y = gap)) +
        geom_errorbar(aes(ymin = gap - sd, ymax = gap + sd), width = 0.18, colour = COL_CURVE, alpha = 0.28, linewidth = 0.45) +
        geom_line(colour = COL_CURVE, linewidth = 1.0) +
        geom_point(colour = COL_CURVE, size = 2.6) +
        geom_point(data = r$gap_df[r$gap_selected_k, ], aes(x = k, y = gap), colour = COL_SELECTED, size = 4.2) +
        scale_x_continuous(breaks = 1:K_MAX) +
        labs(x = "Discrete clusters (K)", y = "Gap statistic", title = sprintf("%s: gap statistic", feature_set_label)) +
        theme_publication()
    }
    p
  }

  # 2x2 grid, not a single column of 4 (too tall/long once rendered at
  # a legible font size) or a single row of 4 (too wide once shrunk to
  # fit a word-processor page). A,B (primary feature set) on top; C,D
  # (sensitivity feature set) on the bottom -- matches the natural
  # primary/sensitivity grouping. Inline grid composer (not theme.R's
  # compose_panels()) -- this script does not source theme.R, matching
  # this project's established R convention that every script is
  # self-contained rather than cross-sourcing another script's file.
  primary_label <- sprintf("Primary feature set (n = %d)", primary_result$n)
  sensitivity_label <- sprintf("Sensitivity feature set (n = %d)", sensitivity_result$n)
  panels <- list(
    list(plot = panel_plot(primary_result, "bic", primary_label), letter = "A"),
    list(plot = panel_plot(primary_result, "gap", primary_label), letter = "B"),
    list(plot = panel_plot(sensitivity_result, "bic", sensitivity_label), letter = "C"),
    list(plot = panel_plot(sensitivity_result, "gap", sensitivity_label), letter = "D")
  )

  panel_letter <- function(letter) {
    grid::textGrob(letter, x = unit(3, "pt"), y = unit(1, "npc") - unit(3, "pt"),
      just = c("left", "top"), gp = grid::gpar(fontface = "bold", fontsize = 28, fontfamily = "Liberation Sans"))
  }

  cairo_pdf(file.path(output_dir_reports, "structure_test.pdf"), width = 15.5, height = 13.0, family = "Liberation Sans")
  grid::grid.newpage()
  layout <- grid::grid.layout(nrow = 2, ncol = 2)
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(panels)) {
    row <- ceiling(i / 2)
    col <- ((i - 1) %% 2) + 1
    grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
    print(panels[[i]]$plot, vp = grid::viewport())
    grid::grid.draw(panel_letter(panels[[i]]$letter))
    grid::popViewport()
  }
  grid::popViewport()
  dev.off()

  cat(sprintf(
    "[OK]   Structure test complete. Final decision: %s. Wrote %s, %s\n",
    final_decision,
    file.path(output_dir_data, "clone_structure_test_results.parquet"),
    file.path(output_dir_reports, "structure_test.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
