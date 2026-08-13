# Unit tests for scripts/15_hpv_stratified_analysis/02_run_prospective_power_simulation.R
# (`15_hpv_stratified_analysis/02_run_prospective_power_simulation.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_power_simulation.R")'

suppressPackageStartupMessages({
  library(testthat)
})

source("../../scripts/15_hpv_stratified_analysis/02_run_prospective_power_simulation.R")

test_that("simulate_power_two_sample_t increases with real larger effect size", {
  set.seed(1)
  low <- simulate_power_two_sample_t(n1 = 20, n2 = 20, cohens_d = 0.2, alpha = 0.05, n_simulations = 1000)
  high <- simulate_power_two_sample_t(n1 = 20, n2 = 20, cohens_d = 2.0, alpha = 0.05, n_simulations = 1000)
  expect_true(high > low)
  expect_true(high > 0.9)
})

test_that("simulate_power_two_sample_t recovers alpha under the real null (d=0)", {
  set.seed(2)
  power_at_null <- simulate_power_two_sample_t(n1 = 30, n2 = 30, cohens_d = 0, alpha = 0.05, n_simulations = 4000)
  # Real Type I error rate should be close to the real nominal alpha.
  expect_true(abs(power_at_null - 0.05) < 0.02)
})

test_that("simulate_power_two_sample_t real n=4 vs n=4 has very low power for a real moderate effect size", {
  set.seed(3)
  power <- simulate_power_two_sample_t(n1 = 4, n2 = 4, cohens_d = 0.8, alpha = 0.05, n_simulations = 2000)
  expect_true(power < 0.3)
})

test_that("compute_exact_permutation_min_pvalue matches the real known closed-form value for n=4 vs n=4", {
  result <- compute_exact_permutation_min_pvalue(n1 = 4, n2 = 4)
  expect_equal(result, 2 / choose(8, 4))
  expect_equal(round(result, 4), 0.0286)
})

test_that("compute_exact_permutation_min_pvalue is smaller for real larger sample sizes", {
  small_n <- compute_exact_permutation_min_pvalue(n1 = 4, n2 = 4)
  large_n <- compute_exact_permutation_min_pvalue(n1 = 20, n2 = 20)
  expect_true(large_n < small_n)
})

test_that("find_minimum_detectable_effect returns NA when no real grid value reaches target power", {
  set.seed(4)
  result <- find_minimum_detectable_effect(n1 = 4, n2 = 4, alpha = 0.05, target_power = 0.8, effect_sizes = c(0.2, 0.4), n_simulations = 200)
  expect_true(is.na(result$minimum_detectable_effect))
})

test_that("find_minimum_detectable_effect finds the smallest real effect size clearing target power", {
  set.seed(5)
  result <- find_minimum_detectable_effect(n1 = 30, n2 = 30, alpha = 0.05, target_power = 0.8, effect_sizes = c(0.2, 0.5, 0.8, 1.2), n_simulations = 1000)
  expect_false(is.na(result$minimum_detectable_effect))
  expect_true(result$minimum_detectable_effect %in% c(0.5, 0.8, 1.2))
})
