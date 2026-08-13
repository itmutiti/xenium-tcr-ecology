# Unit tests for scripts/13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R
# (`13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_variance_partition.R")'

suppressPackageStartupMessages({
  library(testthat)
  library(lme4)
})

source("../../scripts/13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R")

#' Real synthetic data generator with KNOWN, engineered variance
#' components, matching this script's own real nested structure.
make_synthetic_data <- function(patient_sd, clone_sd, residual_sd, n_patients = 8, n_clones_per_patient = 5, n_sections_per_clone = 2, seed = 1) {
  set.seed(seed)
  patient_effect <- rnorm(n_patients, sd = patient_sd)
  rows <- list()
  idx <- 1
  for (p in seq_len(n_patients)) {
    for (c in seq_len(n_clones_per_patient)) {
      clone_effect <- rnorm(1, sd = clone_sd)
      for (s in seq_len(n_sections_per_clone)) {
        rows[[idx]] <- data.frame(
          patient_id = paste0("P", p),
          clone_id = paste0("P", p, "_C", c),
          ecological_structure_score = patient_effect[p] + clone_effect + rnorm(1, sd = residual_sd)
        )
        idx <- idx + 1
      }
    }
  }
  do.call(rbind, rows)
}

test_that("extract_variance_partition pulls all three named real components", {
  data <- make_synthetic_data(patient_sd = 1, clone_sd = 1, residual_sd = 1, seed = 2)
  model <- fit_variance_partition_model(data)
  vc <- extract_variance_partition(model)
  expect_true(vc$patient_var >= 0)
  expect_true(vc$identity_var >= 0)
  expect_true(vc$context_var >= 0)
})

test_that("dominant patient variance is correctly recovered as the largest real component", {
  # Patient effect is engineered to be much larger than clone/residual.
  data <- make_synthetic_data(patient_sd = 5, clone_sd = 0.2, residual_sd = 0.2, seed = 3)
  model <- fit_variance_partition_model(data)
  vc <- extract_variance_partition(model)
  expect_true(vc$patient_var > vc$identity_var)
  expect_true(vc$patient_var > vc$context_var)
})

test_that("dominant clone (identity) variance is correctly recovered as the largest real component", {
  data <- make_synthetic_data(patient_sd = 0.1, clone_sd = 5, residual_sd = 0.1, seed = 4)
  model <- fit_variance_partition_model(data)
  vc <- extract_variance_partition(model)
  expect_true(vc$identity_var > vc$patient_var)
  expect_true(vc$identity_var > vc$context_var)
})

test_that("dominant residual (context) variance is correctly recovered as the largest real component", {
  data <- make_synthetic_data(patient_sd = 0.1, clone_sd = 0.1, residual_sd = 5, seed = 5)
  model <- fit_variance_partition_model(data)
  vc <- extract_variance_partition(model)
  expect_true(vc$context_var > vc$patient_var)
  expect_true(vc$context_var > vc$identity_var)
})

test_that("compute_variance_proportions sums to real 1.0 and matches known ratios", {
  result <- compute_variance_proportions(patient_var = 2, identity_var = 3, context_var = 5)
  expect_equal(result$patient_proportion + result$identity_proportion + result$context_proportion, 1.0)
  expect_equal(result$patient_proportion, 0.2)
  expect_equal(result$identity_proportion, 0.3)
  expect_equal(result$context_proportion, 0.5)
})
