#!/usr/bin/env Rscript
# `14_spatial_interactions_and_barriers/06_benchmark_against_published_barrier_studies.R` -- 06_benchmark_against_published_barrier_studies.R
#
# The registered `q3_literature_benchmark` analysis (governance/analysis_
# registry.tsv, supporting Q3, not independently multiplicity-corrected,
# "a comparison, not a directional hypothesis to be 'confirmed'"): checks
# whether this project's `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R` barrier-topology effect is
# directionally and, where a comparable statistic exists, quantitatively
# concordant with two independently published spatial barrier studies.
#
# Citations, checked against the source publications before use:
#   - Grout JA, Sirven P, Leader AM, et al. "Spatial Positioning and
#     Matrix Programs of Cancer-Associated Fibroblasts Promote T-cell
#     Exclusion in Human Lung Tumors." Cancer Discovery 2022;12(11):
#     2606-2625. DOI 10.1158/2159-8290.CD-21-1714. PMID 36027053.
#     Senior/corresponding authors: Helene Salmon, Ephraim Kenigsberg.
#     Discrepancy from the registry's informal "Grout/Kim 2022" label:
#     no author named "Kim" appears anywhere in this paper's author list
#     (verified via PubMed/PMC directly) -- recorded here rather than
#     silently corrected without comment; this script uses the verified
#     citation.
#   - Hwang WL, Jagadeesh KA, Guo JA, et al. "Single-nucleus and spatial
#     transcriptome profiling of pancreatic cancer identifies
#     multicellular dynamics associated with neoadjuvant treatment."
#     Nature Genetics 2022;54(8):1178-1191. DOI 10.1038/s41588-022-01134-8.
#     PMID 35902743.
#
# Published statistics used for comparison (quoted directly from each
# paper's results, verified via PMC full text before use, not
# fabricated):
#   - Grout et al., Fig 6B: Pearson correlation between stromal alphaSMA
#     (CAF/matrix barrier) density and T-cell infiltration, R = -0.48,
#     p = 1e-10 -- a direct, raw (unadjusted) bivariate analog to this
#     project's barrier-vs-engagement association.
#   - Hwang et al.: malignant-programme hazard ratios for time-to-
#     progression (neural-like progenitor HR=1.62, 95%CI 1.08-2.42,
#     p=0.02; classical HR=0.61, 95%CI 0.46-0.82, p<0.001) and
#     qualitative co-localisation of CD8+ T cells with specific
#     programme-high regions -- a different outcome type (programme-
#     driven prognosis, not barrier-vs-engagement), so this comparison
#     is treated qualitatively only (both studies find that spatially-
#     organised cellular/molecular state non-trivially shapes immune
#     context), per the registry's stated allowance ("Qualitative and,
#     where the published studies report comparable effect-size
#     statistics, quantitative concordance").
#
# Primary output: reports/interactions/literature_benchmark.pdf

suppressPackageStartupMessages({
  library(arrow)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

GROUT_2022_R <- -0.48
GROUT_2022_PVALUE <- 1e-10
GROUT_2022_CITATION <- "Grout et al. 2022, Cancer Discovery 12(11):2606-2625, DOI 10.1158/2159-8290.CD-21-1714"

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

#' Pure, testable: Pearson correlation and p-value between two numeric
#' vectors.
compute_raw_correlation <- function(x, y) {
  result <- stats::cor.test(x, y, method = "pearson")
  list(r = as.numeric(result$estimate), pvalue = result$p.value)
}

#' Pure, testable: direction concordance (same sign) and a magnitude
#' ratio (|project_r| / |published_r|) between this project's
#' correlation and a published one. `magnitude_ratio` is `NA` (not a
#' divide-by-zero artefact) when `published_r == 0`.
classify_concordance <- function(project_r, published_r) {
  list(
    direction_concordant = sign(project_r) == sign(published_r),
    magnitude_ratio = if (published_r == 0) NA_real_ else abs(project_r) / abs(published_r)
  )
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  engagement_path <- file.path(project_root, "data", "derived", "clone_tumour_engagement.parquet")
  barrier_path <- file.path(project_root, "data", "derived", "clone_barrier_metrics.parquet")

  for (p in c(engagement_path, barrier_path)) {
    if (!file.exists(p)) stop(sprintf("'%s' not found.", p))
  }

  engagement <- as.data.frame(read_parquet(engagement_path))
  barrier <- as.data.frame(read_parquet(barrier_path))
  data <- merge(engagement[, c("clone_id", "section_id", "engagement_ratio")], barrier[, c("clone_id", "section_id", "fibroblast_barrier_fraction", "suppressive_myeloid_barrier_fraction")], by = c("clone_id", "section_id"))
  data <- data[stats::complete.cases(data), ]

  cat(sprintf("[INFO] Data: n=%d clone-section rows with a defined barrier fraction.\n", nrow(data)))

  fibroblast_corr <- compute_raw_correlation(data$fibroblast_barrier_fraction, data$engagement_ratio)
  myeloid_corr <- compute_raw_correlation(data$suppressive_myeloid_barrier_fraction, data$engagement_ratio)
  myeloid_concordance <- classify_concordance(myeloid_corr$r, GROUT_2022_R)

  result <- data.frame(
    covariate = c("fibroblast_barrier_fraction", "suppressive_myeloid_barrier_fraction"),
    project_raw_r = c(fibroblast_corr$r, myeloid_corr$r),
    project_raw_pvalue = c(fibroblast_corr$pvalue, myeloid_corr$pvalue),
    published_r = c(NA_real_, GROUT_2022_R),
    published_pvalue = c(NA_real_, GROUT_2022_PVALUE),
    published_citation = c(NA_character_, GROUT_2022_CITATION),
    direction_concordant = c(NA, myeloid_concordance$direction_concordant),
    magnitude_ratio = c(NA_real_, myeloid_concordance$magnitude_ratio)
  )

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "literature_benchmark_results.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "interactions")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # Narrower, taller canvas: a word processor auto-scales a pasted
  # image to fit the page's text width (~6.5in), shrinking every
  # dimension including font size by that ratio. A narrower canvas
  # shrinks less, so the same in-image point size reads larger once
  # pasted.
  open_publication_pdf(file.path(output_dir_reports, "literature_benchmark.pdf"), width = 7.8, height = 7.2)

  plot_data <- data.frame(
    source = factor(c("This study\n(HNSCC, n = 152)", "Grout et al. 2022\n(Fig. 6B)"), levels = c("This study\n(HNSCC, n = 152)", "Grout et al. 2022\n(Fig. 6B)")),
    r = c(myeloid_corr$r, GROUT_2022_R)
  )
  p1 <- ggplot(plot_data, aes(x = source, y = r, fill = source)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = c(PUB_COLORS$primary_analysis, PUB_COLORS$not_significant)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    labs(
      subtitle = sprintf(
        "This study: r = %.3f, p = %.3f\nGrout et al. 2022: r = %.2f, p = %.0e",
        myeloid_corr$r, myeloid_corr$pvalue, GROUT_2022_R, GROUT_2022_PVALUE
      ),
      y = "Raw Pearson r\n(barrier density vs. engagement)", x = NULL
    ) +
    theme_publication()
  print(p1)

  dev.off()

  cat(sprintf("[INFO] Raw correlations: fibroblast r=%.4f (p=%.4f), suppressive_myeloid r=%.4f (p=%.4f)\n", fibroblast_corr$r, fibroblast_corr$pvalue, myeloid_corr$r, myeloid_corr$pvalue))
  cat(sprintf(
    "[INFO] vs. %s: r=%.2f (p=%.0e) -- direction %s, magnitude ratio (this project / published) = %.3f\n",
    GROUT_2022_CITATION, GROUT_2022_R, GROUT_2022_PVALUE,
    if (myeloid_concordance$direction_concordant) "concordant" else "discordant", myeloid_concordance$magnitude_ratio
  ))
  cat(sprintf(
    "[OK]   Literature benchmark complete. Wrote %s, %s\n",
    file.path(output_dir_data, "literature_benchmark_results.parquet"),
    file.path(output_dir_reports, "literature_benchmark.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
