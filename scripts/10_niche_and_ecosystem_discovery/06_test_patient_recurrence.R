#!/usr/bin/env Rscript
# `10_niche_and_ecosystem_discovery/06_test_patient_recurrence.R` -- 06_test_patient_recurrence.R
#
# Tests whether each `10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py` ecosystem's abundance is a recurrent,
# patient-level property -- consistent across a patient's technical
# replicate Xenium runs -- or dominated by replicate-to-replicate
# technical noise, via a random-intercept mixed model with
# technical-replicate sections nested within patient
# (metadata/sample_manifest.tsv: 7 of 11 patients have 2 replicate runs,
# 4 have 1; `is_technical_replicate` / `patient_id` are existing
# columns, not derived by parsing `section_id` strings).
#
# "Recurrence" here means statistical/spatial recurrence of ecosystem
# structure across the cohort -- matching `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s "recurrent
# neighbourhoods" terminology -- not the clinical `recurrence_status`
# (Recurrent/Primary) field in the sample manifest, which is an unrelated
# clinical outcome and explicitly out of scope here, same as this
# script's "HPV is not tested here" scope note: no clinical-outcome
# association is tested in this milestone either.
#
# Primary output: reports/niches/recurrence_models.pdf

suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(lme4)
  library(ggplot2)
})

PROJECT_ROOT_MARKER <- "manifests/project_paths.yaml"
ENV_ROOT_VAR <- "XENIUM_TCR_ECOLOGY_ROOT"

N_BOOTSTRAP <- 500

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

#' Fits the recurrence model for one ecosystem's per-section abundance
#' data: a random-intercept mixed model of logit(abundance) on patient
#' identity, with replicate sections nested within patient.
fit_recurrence_model <- function(data) {
  lme4::lmer(logit_abundance ~ 1 + (1 | patient_id), data = data, REML = TRUE)
}

#' Pure, testable: extracts the patient-level and residual variance
#' components from a fitted `fit_recurrence_model` object.
extract_variance_components <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  list(
    patient_var = vc$vcov[vc$grp == "patient_id"],
    residual_var = vc$vcov[vc$grp == "Residual"]
  )
}

#' Pure, testable: intraclass correlation coefficient (ICC) -- the
#' fraction of total variance attributable to patient identity. ICC near
#' 1 means an ecosystem's abundance is a strongly recurrent, patient-
#' specific property (technical replicates agree closely); ICC near 0
#' means abundance is dominated by replicate-to-replicate technical
#' noise, not patient-level structure.
compute_icc <- function(patient_var, residual_var) {
  patient_var / (patient_var + residual_var)
}

main <- function() {
  project_root <- find_project_root(parse_project_root_arg())
  source(file.path(project_root, "src", "xenium_tcr_ecology", "viz", "theme.R"))

  metrics_path <- file.path(project_root, "data", "derived", "ecosystem_metrics.parquet")
  manifest_path <- file.path(project_root, "metadata", "sample_manifest.tsv")
  if (!file.exists(metrics_path)) stop(sprintf("'%s' not found. Run `10_niche_and_ecosystem_discovery/05_quantify_ecosystem_abundance_and_topology.py` first.", metrics_path))
  if (!file.exists(manifest_path)) stop(sprintf("'%s' not found.", manifest_path))

  metrics <- as.data.frame(read_parquet(metrics_path))
  manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE)

  merged <- merge(metrics, manifest[, c("section_id", "patient_id")], by = "section_id")
  if (nrow(merged) != nrow(metrics)) {
    stop(sprintf("Merge with sample manifest changed row count (%d -> %d) -- check for unmatched section_id values.", nrow(metrics), nrow(merged)))
  }
  merged$logit_abundance <- qlogis(merged$abundance)

  ecosystem_labels <- sort(unique(merged$ecosystem_label))
  results <- vector("list", length(ecosystem_labels))
  models <- list()

  for (i in seq_along(ecosystem_labels)) {
    label <- ecosystem_labels[i]
    sub <- merged[merged$ecosystem_label == label, ]
    model <- fit_recurrence_model(sub)
    models[[label]] <- model
    vc <- extract_variance_components(model)
    icc <- compute_icc(vc$patient_var, vc$residual_var)

    set.seed(RNG_SEED)
    boot <- tryCatch(
      lme4::bootMer(
        model,
        FUN = function(m) {
          vc_b <- extract_variance_components(m)
          compute_icc(vc_b$patient_var, vc_b$residual_var)
        },
        nsim = N_BOOTSTRAP, use.u = FALSE, type = "parametric"
      ),
      error = function(e) NULL
    )
    if (!is.null(boot)) {
      ci <- stats::quantile(boot$t, c(0.025, 0.975), na.rm = TRUE)
    } else {
      ci <- c(NA_real_, NA_real_)
    }

    results[[i]] <- data.frame(
      ecosystem_label = label,
      n_sections = nrow(sub),
      n_patients = length(unique(sub$patient_id)),
      patient_var = vc$patient_var,
      residual_var = vc$residual_var,
      icc = icc,
      icc_ci_low = unname(ci[1]),
      icc_ci_high = unname(ci[2])
    )
  }
  result <- do.call(rbind, results)
  rownames(result) <- NULL

  output_dir_data <- file.path(project_root, "data", "derived")
  write_parquet(result, file.path(output_dir_data, "ecosystem_recurrence_models.parquet"))

  output_dir_reports <- file.path(project_root, "reports", "niches")
  dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)
  # A on top, B (a 6-ecosystem facet grid) below -- but B's facets run
  # 3-per-row (ncol=3, 2 rows) rather than 2-per-row (3 rows), keeping
  # this pairing close to square instead of the excessively tall canvas
  # a taller facet grid forces. A word processor auto-scales a pasted
  # image to fit the page's text width (~6.5in), shrinking every
  # dimension including font size by that ratio, so an oversized height
  # is pure wasted legibility, not just an aesthetic issue.
  open_publication_pdf(file.path(output_dir_reports, "recurrence_models.pdf"), width = 18.0, height = 13.0)

  result_ordered <- result[order(result$icc), ]
  result_ordered$ecosystem_label <- factor(result_ordered$ecosystem_label, levels = result_ordered$ecosystem_label)
  p1 <- ggplot(result_ordered, aes(x = icc, y = ecosystem_label)) +
    geom_pointrange(aes(xmin = icc_ci_low, xmax = icc_ci_high), colour = PUB_COLORS$primary_analysis, linewidth = 1.0, size = 0.65) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = PUB_COLORS$reference_line, linewidth = 0.7) +
    xlim(0, 1) +
    labs(
      subtitle = sprintf("Patient-level ICC, 95%% CI\nn = %d sections, %d patients", length(unique(metrics$section_id)), length(unique(merged$patient_id))),
      x = "ICC (variance attributable to patient identity)", y = NULL
    ) +
    theme_publication()

  merged_ordered <- merge(merged, manifest[, c("section_id", "run_number")], by = "section_id")
  p2 <- ggplot(merged_ordered, aes(x = reorder(patient_id, abundance, FUN = median), y = abundance, colour = factor(run_number))) +
    geom_point(size = 2.4, alpha = 0.85) +
    facet_wrap(~ecosystem_label, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c(OKABE_ITO$blue, OKABE_ITO$vermillion, OKABE_ITO$bluish_green)) +
    labs(
      subtitle = "Ecosystem abundance by patient, coloured by technical replicate run",
      x = "Patient (ordered by median abundance)", y = "Abundance", colour = "Run"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
      panel.spacing.y = unit(1.4, "lines"), panel.spacing.x = unit(1.2, "lines")
    )

  compose_panels(list(p1, p2), ncol = 1, heights = c(1, 1.6))
  dev.off()

  cat("[INFO] ICC (patient-level recurrence) by ecosystem:\n")
  for (i in seq_len(nrow(result))) {
    cat(sprintf(
      "[INFO]   %-32s ICC=%.3f (95%% CI %.3f-%.3f)\n",
      result$ecosystem_label[i], result$icc[i], result$icc_ci_low[i], result$icc_ci_high[i]
    ))
  }
  cat(sprintf(
    "[OK]   %d ecosystem(s) modelled across %d sections, %d patients. Wrote %s, %s\n",
    length(ecosystem_labels), length(unique(metrics$section_id)), length(unique(merged$patient_id)),
    file.path(output_dir_data, "ecosystem_recurrence_models.parquet"),
    file.path(output_dir_reports, "recurrence_models.pdf")
  ))
}

if (sys.nframe() == 0) {
  main()
}
