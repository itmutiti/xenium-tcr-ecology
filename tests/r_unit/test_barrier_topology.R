# Unit tests for scripts/14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R
# (`14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`, prespecified q3_barrier_topology_confirmatory). Run from
# the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_barrier_topology.R")'

suppressPackageStartupMessages({
  library(testthat)
  library(lme4)
})

source("../../scripts/14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R")
source("../../scripts/14_spatial_interactions_and_barriers/07_ablate_covariates_for_barrier_effect.R")

#' Real synthetic (patient_id, clone_id) key generator with 2 real
#' replicate rows per clone (matching this project's own real
#' clone-section structure, e.g. n=261 rows / 158 clones) -- a single
#' row per clone makes the clone-level random effect non-identifiable
#' (`lme4` correctly refuses), the same real bug already caught and
#' fixed in `13_clone_ecology_confirmatory_models/03_test_clone_size_as_confounder.R`'s own synthetic test data.
make_keys <- function(n_patients, n_clones_per_patient, n_reps_per_clone = 2) {
  patient_id <- rep(paste0("P", seq_len(n_patients)), each = n_clones_per_patient * n_reps_per_clone)
  clone_id <- unlist(lapply(seq_len(n_patients), function(p) rep(paste0("P", p, "_C", seq_len(n_clones_per_patient)), each = n_reps_per_clone)))
  data.frame(patient_id = patient_id, clone_id = clone_id)
}

test_that("compute_niche_composition_wide gives real per-(clone,section) archetype fractions summing to 1", {
  domain_data <- data.frame(
    clone_id = c("C1", "C1", "C1", "C1", "C2", "C2"),
    section_id = c("S1", "S1", "S1", "S1", "S1", "S1"),
    archetype = c(1, 1, 2, 3, 2, 2)
  )
  result <- compute_niche_composition_wide(domain_data, levels = 1:3)
  c1 <- result[result$clone_id == "C1" & result$section_id == "S1", ]
  expect_equal(c1$niche_archetype_1_fraction, 0.5)
  expect_equal(c1$niche_archetype_2_fraction, 0.25)
  expect_equal(c1$niche_archetype_3_fraction, 0.25)
  c2 <- result[result$clone_id == "C2" & result$section_id == "S1", ]
  expect_equal(c2$niche_archetype_2_fraction, 1.0)
  expect_equal(c2$niche_archetype_1_fraction, 0.0)
  expect_equal(rowSums(result[, c("niche_archetype_1_fraction", "niche_archetype_2_fraction", "niche_archetype_3_fraction")]), c(1, 1))
})

test_that("compute_marginal_r2 recovers a known real value on an engineered model", {
  set.seed(10)
  keys <- make_keys(n_patients = 10, n_clones_per_patient = 10)
  n <- nrow(keys)
  x <- rnorm(n)
  clone_effect <- rnorm(length(unique(keys$clone_id)), sd = 0.05)[match(keys$clone_id, unique(keys$clone_id))]
  patient_effect <- rnorm(length(unique(keys$patient_id)), sd = 0.05)[match(keys$patient_id, unique(keys$patient_id))]
  # Fixed-effect signal dominates; random/residual noise is small real
  # relative to it -- marginal R2 should be close to 1.
  y <- 5 * x + rnorm(n, sd = 0.05) + clone_effect + patient_effect
  data <- data.frame(keys, x = x, y = y)
  model <- lmer(y ~ x + (1 | patient_id / clone_id), data = data, REML = FALSE)
  r2 <- compute_marginal_r2(model)
  expect_true(r2 > 0.9)
  expect_true(r2 <= 1.0)
})

test_that("compute_marginal_r2 is low when the fixed effect carries real no signal", {
  set.seed(11)
  keys <- make_keys(n_patients = 10, n_clones_per_patient = 10)
  n <- nrow(keys)
  x <- rnorm(n)
  clone_effect <- rnorm(length(unique(keys$clone_id)), sd = 1)[match(keys$clone_id, unique(keys$clone_id))]
  patient_effect <- rnorm(length(unique(keys$patient_id)), sd = 1)[match(keys$patient_id, unique(keys$patient_id))]
  y <- rnorm(n, sd = 1) + clone_effect + patient_effect
  data <- data.frame(keys, x = x, y = y)
  model <- lmer(y ~ x + (1 | patient_id / clone_id), data = data, REML = FALSE)
  r2 <- compute_marginal_r2(model)
  expect_true(r2 < 0.1)
})

test_that("compare_nested_models recovers a significant LRT when the added term carries real engineered signal", {
  set.seed(12)
  keys <- make_keys(n_patients = 12, n_clones_per_patient = 8)
  n <- nrow(keys)
  covariate <- rnorm(n)
  barrier <- rnorm(n)
  clone_effect <- rnorm(length(unique(keys$clone_id)), sd = 0.3)[match(keys$clone_id, unique(keys$clone_id))]
  patient_effect <- rnorm(length(unique(keys$patient_id)), sd = 0.3)[match(keys$patient_id, unique(keys$patient_id))]
  y <- 0.2 * covariate + 3 * barrier + rnorm(n, sd = 0.3) + clone_effect + patient_effect
  data <- data.frame(keys, covariate = covariate, barrier = barrier, y = y)
  baseline <- lmer(y ~ covariate + (1 | patient_id / clone_id), data = data, REML = FALSE)
  full <- lmer(y ~ covariate + barrier + (1 | patient_id / clone_id), data = data, REML = FALSE)
  comparison <- compare_nested_models(baseline, full)
  expect_true(comparison$pvalue < 0.05)
  expect_true(comparison$chisq > 0)
})

test_that("compare_nested_models fails to reject when the added term carries real no signal", {
  set.seed(13)
  keys <- make_keys(n_patients = 12, n_clones_per_patient = 8)
  n <- nrow(keys)
  covariate <- rnorm(n)
  barrier <- rnorm(n)
  clone_effect <- rnorm(length(unique(keys$clone_id)), sd = 1)[match(keys$clone_id, unique(keys$clone_id))]
  patient_effect <- rnorm(length(unique(keys$patient_id)), sd = 1)[match(keys$patient_id, unique(keys$patient_id))]
  y <- 0.2 * covariate + rnorm(n, sd = 1) + clone_effect + patient_effect
  data <- data.frame(keys, covariate = covariate, barrier = barrier, y = y)
  baseline <- lmer(y ~ covariate + (1 | patient_id / clone_id), data = data, REML = FALSE)
  full <- lmer(y ~ covariate + barrier + (1 | patient_id / clone_id), data = data, REML = FALSE)
  comparison <- compare_nested_models(baseline, full)
  expect_true(comparison$pvalue > 0.05)
})

test_that("fit_and_extract_barrier_effect recovers a known engineered barrier coefficient with no adjustment covariates", {
  set.seed(14)
  keys <- make_keys(n_patients = 12, n_clones_per_patient = 8)
  n <- nrow(keys)
  fibroblast_barrier_fraction <- runif(n)
  suppressive_myeloid_barrier_fraction <- runif(n)
  clone_effect <- rnorm(length(unique(keys$clone_id)), sd = 0.05)[match(keys$clone_id, unique(keys$clone_id))]
  patient_effect <- rnorm(length(unique(keys$patient_id)), sd = 0.05)[match(keys$patient_id, unique(keys$patient_id))]
  engagement_ratio <- -2 * suppressive_myeloid_barrier_fraction + rnorm(n, sd = 0.05) + clone_effect + patient_effect
  data <- data.frame(keys, fibroblast_barrier_fraction, suppressive_myeloid_barrier_fraction, engagement_ratio)
  result <- fit_and_extract_barrier_effect(data, character(0), "barrier_only")
  expect_equal(result$step, "barrier_only")
  expect_equal(result$n_adjustment_covariates, 0)
  expect_true(result$estimate < -1)  # real, strong engineered negative effect recovered
  expect_true(result$p_value < 0.05)
  expect_true(result$ci_low < result$estimate && result$ci_high > result$estimate)
})

test_that("fit_and_extract_barrier_effect gives a real non-significant estimate when the barrier carries no engineered signal", {
  set.seed(15)
  keys <- make_keys(n_patients = 12, n_clones_per_patient = 8)
  n <- nrow(keys)
  fibroblast_barrier_fraction <- runif(n)
  suppressive_myeloid_barrier_fraction <- runif(n)
  covariate_x <- rnorm(n)
  clone_effect <- rnorm(length(unique(keys$clone_id)), sd = 1)[match(keys$clone_id, unique(keys$clone_id))]
  patient_effect <- rnorm(length(unique(keys$patient_id)), sd = 1)[match(keys$patient_id, unique(keys$patient_id))]
  engagement_ratio <- 0.1 * covariate_x + rnorm(n, sd = 1) + clone_effect + patient_effect
  data <- data.frame(keys, fibroblast_barrier_fraction, suppressive_myeloid_barrier_fraction, covariate_x, engagement_ratio)
  result <- fit_and_extract_barrier_effect(data, "covariate_x", "+ covariate_x")
  expect_equal(result$n_adjustment_covariates, 1)
  expect_true(result$p_value > 0.05)
})
