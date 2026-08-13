#!/usr/bin/env Rscript
# `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R` -- 02_discover_neighbourhood_archetypes.R
#
# Identifies recurrent neighbourhoods ("archetypes") from `10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`'s
# local composition vectors via consensus clustering (Monti et al. 2003),
# and quantifies clustering stability across candidate K via the
# Proportion of Ambiguous Clustering (PAC) statistic (Senbabaoglu et al.
# 2014) -- lower PAC means a more decisive, reproducible partition.
#
# Uses only the radius_30.0um composition columns (not all three scales
# pooled) -- `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py` calibrated 30um as this project's primary spatial
# scale (config/graph_parameters.yaml: calibrated_radius_um: 30.0), so
# archetype discovery reuses that same calibrated scale rather than
# introducing a new, independently-chosen multi-scale feature space.
#
# Consensus clustering runs on a fixed-size data subsample (stability
# analysis is O(n^2) in the consensus matrix; the final archetype
# assignment is fit on the full dataset once the stable K is selected).
#
# Primary output: data/derived/neighbourhood_archetypes.parquet
# Report: reports/niches/archetype_stability.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

CANDIDATE_K <- c(4, 6, 8, 10, 12) # spans coarser-than- to as-fine-as `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s 12 major lineages
N_STABILITY_CELLS <- 3000 # consensus matrix is O(n^2); 3000^2 = 9e6 entries, feasible in base R
N_BOOTSTRAP <- 50
BOOTSTRAP_SUBSAMPLE_FRAC <- 0.8
PAC_LOWER <- 0.1
PAC_UPPER <- 0.9
FINAL_KMEANS_NSTART <- 25

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

#' Consensus (co-clustering) matrix via bootstrap resampled k-means
#' (Monti et al. 2003). Pure, testable: takes a plain numeric matrix.
#' Returns an n x n matrix where entry [i,j] is the fraction of bootstrap
#' resamples containing BOTH i and j that placed them in the same cluster
#' (NA where i and j were never co-sampled).
compute_consensus_matrix <- function(data, k, n_bootstrap = N_BOOTSTRAP, subsample_frac = BOOTSTRAP_SUBSAMPLE_FRAC, seed = RNG_SEED) {
  n <- nrow(data)
  same_cluster_count <- matrix(0L, n, n)
  co_sampled_count <- matrix(0L, n, n)

  for (b in seq_len(n_bootstrap)) {
    set.seed(seed + b)
    idx <- sort(sample(n, size = round(n * subsample_frac)))
    km <- kmeans(data[idx, , drop = FALSE], centers = k, nstart = 5, iter.max = 100)
    same <- outer(km$cluster, km$cluster, `==`)
    same_cluster_count[idx, idx] <- same_cluster_count[idx, idx] + same
    co_sampled_count[idx, idx] <- co_sampled_count[idx, idx] + 1L
  }

  consensus <- same_cluster_count / co_sampled_count
  diag(consensus) <- NA
  consensus
}

#' Proportion of Ambiguous Clustering (Senbabaoglu et al. 2014): the
#' fraction of consensus-matrix entries falling in the "ambiguous" middle
#' zone [pac_lower, pac_upper], excluding never-co-sampled (NA) pairs.
#' Lower PAC = a more decisive (bimodal, near-0/near-1) consensus matrix,
#' i.e. a more stable clustering at that K.
compute_pac_score <- function(consensus_matrix, pac_lower = PAC_LOWER, pac_upper = PAC_UPPER) {
  values <- consensus_matrix[upper.tri(consensus_matrix)]
  values <- values[!is.na(values)]
  mean(values >= pac_lower & values <= pac_upper)
}

#' Selects the K with the lowest (best) PAC score. `pac_by_k` must be a
#' named numeric vector, names are the candidate K values as strings.
select_stable_k <- function(pac_by_k) {
  as.integer(names(pac_by_k)[which.min(pac_by_k)])
}

#' For each archetype cluster, the lineage column with the highest mean
#' composition fraction -- a direct interpretability label ("this
#' archetype is Fibroblast-dominant"), not an arbitrary cluster index.
label_archetypes_by_dominant_lineage <- function(centroids, lineage_names) {
  apply(centroids[, lineage_names, drop = FALSE], 1, function(row) lineage_names[which.max(row)])
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  compositions_path <- file.path(project_root, "data", "derived", "local_compositions.parquet")
  if (!file.exists(compositions_path)) {
    stop(sprintf("'%s' not found. Run `10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py` first.", compositions_path))
  }
  compositions <- as.data.frame(read_parquet(compositions_path))

  scale_prefix <- "radius_30.0um__"
  composition_cols <- grep(paste0("^", scale_prefix), colnames(compositions), value = TRUE)
  if (length(composition_cols) == 0) {
    stop(sprintf("No columns matching '%s*' found in '%s'.", scale_prefix, compositions_path))
  }
  lineage_names <- sub(paste0("^", scale_prefix), "", composition_cols)

  n_total <- nrow(compositions)
  complete <- stats::complete.cases(compositions[, composition_cols])
  n_zero_degree <- sum(!complete)
  fitted <- compositions[complete, ]
  feature_matrix <- as.matrix(fitted[, composition_cols])
  colnames(feature_matrix) <- lineage_names

  # --- Stage 1: consensus-clustering stability analysis on a subsample ---
  set.seed(RNG_SEED)
  n_stability <- min(N_STABILITY_CELLS, nrow(feature_matrix))
  stability_idx <- sample(nrow(feature_matrix), size = n_stability)
  stability_matrix <- feature_matrix[stability_idx, , drop = FALSE]

  pac_by_k <- setNames(numeric(length(CANDIDATE_K)), as.character(CANDIDATE_K))
  for (k in CANDIDATE_K) {
    consensus <- compute_consensus_matrix(stability_matrix, k = k)
    pac_by_k[as.character(k)] <- compute_pac_score(consensus)
  }
  stable_k <- select_stable_k(pac_by_k)
  cat("[INFO] PAC score by candidate K (lower = more stable):\n")
  for (k in CANDIDATE_K) cat(sprintf("[INFO]   K=%2d  PAC=%.4f\n", k, pac_by_k[as.character(k)]))
  cat(sprintf("[INFO] Selected K=%d (minimum PAC).\n", stable_k))

  # --- Stage 2: final archetype assignment on the full dataset ---
  set.seed(RNG_SEED)
  final_km <- kmeans(feature_matrix, centers = stable_k, nstart = FINAL_KMEANS_NSTART, iter.max = 200)

  centroids <- as.data.frame(final_km$centers)
  centroids$archetype <- seq_len(stable_k)
  centroids$dominant_lineage <- label_archetypes_by_dominant_lineage(final_km$centers, lineage_names)
  centroids$n_cells <- as.integer(table(factor(final_km$cluster, levels = seq_len(stable_k))))

  # `fitted$cell_id` is the actual cell barcode (arrow preserves pandas'
  # named index as a genuine column) -- not `rownames(fitted)`, which is
  # just R's default 1..n row-label sequence and would silently discard
  # the cell identity needed for downstream joins.
  result <- fitted[, c("cell_id", "section_id"), drop = FALSE]
  result$archetype <- final_km$cluster
  result$dominant_lineage <- centroids$dominant_lineage[final_km$cluster]

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "neighbourhood_archetypes.parquet"))
  write_parquet(centroids, file.path(output_dir_data, "neighbourhood_archetype_centroids.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "niches")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Side-by-side, not stacked: for 2 panels a single moderate-width
  # row keeps the aspect ratio closer to square than either a wide
  # 2-in-a-row or a tall 2x1 stack would. B gets more width than A
  # (a 10-lineage heatmap with a colour legend vs. a simple line plot).
  open_publication_pdf(file.path(output_dir_reports, "archetype_stability.pdf"), width = 16.5, height = 7.5)

  pac_df <- data.frame(K = CANDIDATE_K, PAC = as.numeric(pac_by_k))
  p1 <- ggplot(pac_df, aes(x = K, y = PAC)) +
    geom_line(colour = PUB_COLORS$primary_analysis, linewidth = 1.2) +
    geom_point(size = 2.8, colour = PUB_COLORS$primary_analysis) +
    geom_point(data = pac_df[pac_df$K == stable_k, ], aes(x = K, y = PAC), colour = PUB_COLORS$sensitivity_analysis, size = 4.6) +
    scale_x_continuous(breaks = CANDIDATE_K) +
    labs(
      subtitle = sprintf("Consensus clustering stability; selected K = %d", stable_k),
      x = "Number of archetypes (K)", y = "Proportion of ambiguous clustering"
    ) +
    theme_publication()

  centroid_long <- data.frame(
    archetype = rep(centroids$archetype, times = length(lineage_names)),
    lineage = rep(lineage_names, each = nrow(centroids)),
    fraction = as.vector(as.matrix(centroids[, lineage_names]))
  )
  centroid_long$archetype_label <- sprintf("A%d (%s, n = %d)", centroid_long$archetype,
    centroids$dominant_lineage[match(centroid_long$archetype, centroids$archetype)],
    centroids$n_cells[match(centroid_long$archetype, centroids$archetype)])
  p2 <- ggplot(centroid_long, aes(x = lineage, y = archetype_label, fill = fraction)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Mean\nfraction") +
    labs(
      subtitle = sprintf("Archetype centroids (K = %d)\n30 um neighbourhood composition", stable_k),
      x = "Neighbour lineage", y = NULL
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  compose_panels(list(p1, p2), ncol = 2, widths = c(1, 1.5))
  dev.off()

  cat(sprintf(
    "[OK]   %d cells (%d zero-degree at radius_30.0um excluded), K=%d archetypes.\n",
    nrow(fitted), n_zero_degree, stable_k
  ))
  cat(sprintf(
    "[OK]   Wrote %s, %s, %s\n",
    file.path(output_dir_data, "neighbourhood_archetypes.parquet"),
    file.path(output_dir_data, "neighbourhood_archetype_centroids.parquet"),
    file.path(output_dir_reports, "archetype_stability.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
