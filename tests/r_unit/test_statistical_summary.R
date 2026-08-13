# Unit tests for scripts/17_statistical_closure_and_release/01_control_multiplicity_and_report_effects.R
# (`17_statistical_closure_and_release/01_control_multiplicity_and_report_effects.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_statistical_summary.R")'

suppressPackageStartupMessages({
  library(testthat)
})

source("../../scripts/17_statistical_closure_and_release/01_control_multiplicity_and_report_effects.R")

test_that("compute_wald_pvalue recovers a known real value for a real z=2 case", {
  result <- compute_wald_pvalue(estimate = 2.0, se = 1.0)
  expect_equal(round(result, 4), round(2 * (1 - pnorm(2)), 4))
})

test_that("compute_wald_pvalue gives p near 1 for a real null estimate", {
  result <- compute_wald_pvalue(estimate = 0.0, se = 1.0)
  expect_equal(result, 1.0)
})

test_that("apply_bonferroni multiplies and caps at 1.0", {
  result <- apply_bonferroni(c(0.01, 0.5, 0.001), n_family = 5)
  expect_equal(result, c(0.05, 1.0, 0.005))
})

test_that("apply_bonferroni caps a large real product at exactly 1.0", {
  result <- apply_bonferroni(c(0.9), n_family = 5)
  expect_equal(result, 1.0)
})
