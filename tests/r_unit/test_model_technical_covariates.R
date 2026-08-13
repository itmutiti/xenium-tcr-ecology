# Unit tests for scripts/05_preprocessing_and_normalisation/04_model_technical_covariates.R
# (`05_preprocessing_and_normalisation/04_model_technical_covariates.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_model_technical_covariates.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/05_preprocessing_and_normalisation/04_model_technical_covariates.R")

test_that("decompose_variance recovers a known variance structure from simulated data", {
  set.seed(1)
  n_patients <- 8
  n_sections_per_patient <- 2
  n_per_section <- 200

  patient_effect_sd <- 2.0
  section_effect_sd <- 0.5
  residual_sd <- 1.0
  fixed_beta <- 1.5

  rows <- list()
  for (p in seq_len(n_patients)) {
    patient_re <- rnorm(1, sd = patient_effect_sd)
    for (s in seq_len(n_sections_per_patient)) {
      section_re <- rnorm(1, sd = section_effect_sd)
      depth_z <- rnorm(n_per_section)
      y <- patient_re + section_re + fixed_beta * depth_z + rnorm(n_per_section, sd = residual_sd)
      rows[[length(rows) + 1]] <- data.frame(
        patient_id = paste0("P", p),
        section_id = paste0("P", p, "_S", s),
        depth_z = depth_z,
        control_burden_z = rnorm(n_per_section),
        y = y
      )
    }
  }
  df <- do.call(rbind, rows)

  model <- lme4::lmer(y ~ depth_z + control_burden_z + (1 | patient_id) + (1 | section_id), data = df)
  result <- decompose_variance(model)

  # Patient variance (sd=2.0, var=4.0) should be the largest single random-
  # effect component, clearly bigger than section variance (sd=0.5, var=0.25).
  patient_frac <- result$fraction[result$component == "patient"]
  section_frac <- result$fraction[result$component == "section (run nested in patient)"]
  fixed_frac <- result$fraction[result$component == "fixed (depth + control burden)"]

  expect_gt(patient_frac, section_frac)
  expect_gt(patient_frac, 0.3)   # patient variance dominates by construction
  expect_gt(fixed_frac, 0.05)    # fixed effect (beta=1.5, non-trivial) has real signal
  expect_equal(sum(result$fraction), 1.0, tolerance = 1e-6)
})

test_that("decompose_variance fractions are always non-negative and sum to 1", {
  set.seed(2)
  n <- 400
  df <- data.frame(
    patient_id = rep(paste0("P", 1:10), each = n / 10),
    section_id = rep(paste0("S", 1:20), each = n / 20),
    depth_z = rnorm(n),
    control_burden_z = rnorm(n),
    y = rnorm(n)
  )
  model <- lme4::lmer(y ~ depth_z + control_burden_z + (1 | patient_id) + (1 | section_id), data = df)
  result <- decompose_variance(model)

  expect_true(all(result$fraction >= 0))
  expect_equal(sum(result$fraction), 1.0, tolerance = 1e-6)
  expect_equal(nrow(result), 4)
})
