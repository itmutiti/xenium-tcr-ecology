# Unit tests for scripts/10_niche_and_ecosystem_discovery/06_test_patient_recurrence.R
# (`10_niche_and_ecosystem_discovery/06_test_patient_recurrence.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_patient_recurrence.R")'

suppressPackageStartupMessages({
  library(testthat)
  library(lme4)
})

source("../../scripts/10_niche_and_ecosystem_discovery/06_test_patient_recurrence.R")

test_that("compute_icc is the patient-variance share of total variance", {
  expect_equal(compute_icc(patient_var = 4, residual_var = 1), 0.8)
  expect_equal(compute_icc(patient_var = 0, residual_var = 1), 0.0)
  expect_equal(compute_icc(patient_var = 1, residual_var = 1), 0.5)
})

test_that("extract_variance_components pulls patient and residual variance from a fitted model", {
  set.seed(1)
  n_patients <- 10
  patient_id <- rep(paste0("P", seq_len(n_patients)), each = 3)
  patient_effect <- rep(rnorm(n_patients, sd = 2), each = 3)
  y <- patient_effect + rnorm(length(patient_id), sd = 0.1)
  data <- data.frame(logit_abundance = y, patient_id = patient_id)

  model <- fit_recurrence_model(data)
  vc <- extract_variance_components(model)

  expect_true(vc$patient_var > 0)
  expect_true(vc$residual_var > 0)
  # Strong patient signal (sd=2) vs tiny residual noise (sd=0.1) -> patient
  # variance should dominate total variance by a wide margin.
  expect_true(vc$patient_var > vc$residual_var)
})

test_that("high between-patient, low within-patient variance gives ICC near 1", {
  set.seed(2)
  n_patients <- 12
  patient_id <- rep(paste0("P", seq_len(n_patients)), each = 2)
  patient_effect <- rep(rnorm(n_patients, sd = 5), each = 2)
  y <- patient_effect + rnorm(length(patient_id), sd = 0.05)
  data <- data.frame(logit_abundance = y, patient_id = patient_id)

  model <- fit_recurrence_model(data)
  vc <- extract_variance_components(model)
  icc <- compute_icc(vc$patient_var, vc$residual_var)

  expect_true(icc > 0.9)
})

test_that("no between-patient variance, only within-patient noise gives ICC near 0", {
  set.seed(3)
  n_patients <- 12
  patient_id <- rep(paste0("P", seq_len(n_patients)), each = 4)
  y <- rnorm(length(patient_id), sd = 1) # no real patient_effect term at all
  data <- data.frame(logit_abundance = y, patient_id = patient_id)

  model <- fit_recurrence_model(data)
  vc <- extract_variance_components(model)
  icc <- compute_icc(vc$patient_var, vc$residual_var)

  expect_true(icc < 0.3)
})
