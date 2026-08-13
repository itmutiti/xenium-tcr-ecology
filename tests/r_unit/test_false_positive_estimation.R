# Unit tests for scripts/08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R
# (`08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_false_positive_estimation.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R")

test_that("compute_background_rate averages the two independent controls", {
  expect_equal(compute_background_rate(0.02, 0.04), 0.03)
  expect_equal(compute_background_rate(0.0, 0.0), 0.0)
})

test_that("compute_empirical_fpr is the expected-false-positives fraction", {
  # 1000 cells at background rate 0.01 -> 10 expected false positives.
  # 20 observed positives -> fpr = 10/20 = 0.5.
  fpr <- compute_empirical_fpr(own_n_tcells = 1000, own_n_detected = 20, background_rate = 0.01)
  expect_equal(fpr, 0.5)
})

test_that("compute_empirical_fpr is capped at 1.0", {
  # Background alone would predict more false positives than were
  # observed -- must report "up to 100%", not a value above 1.0.
  fpr <- compute_empirical_fpr(own_n_tcells = 1000, own_n_detected = 5, background_rate = 0.01)
  expect_equal(fpr, 1.0)
})

test_that("compute_empirical_fpr is near zero when background rate is near zero", {
  fpr <- compute_empirical_fpr(own_n_tcells = 1000, own_n_detected = 20, background_rate = 0.0)
  expect_equal(fpr, 0.0)
})

test_that("compute_empirical_fpr is vectorised over probes", {
  fpr <- compute_empirical_fpr(
    own_n_tcells = c(1000, 1000),
    own_n_detected = c(20, 5),
    background_rate = c(0.01, 0.01)
  )
  expect_equal(fpr, c(0.5, 1.0))
})
