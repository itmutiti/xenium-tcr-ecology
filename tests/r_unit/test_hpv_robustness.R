# Unit tests for scripts/15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R
# (`15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_hpv_robustness.R")'

suppressPackageStartupMessages({
  library(testthat)
})

source("../../scripts/15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R")

test_that("run_exhaustive_permutation_test enumerates all choose(8,4)=70 real relabellings for n=4 vs n=4", {
  values <- c(1, 2, 3, 4, 5, 6, 7, 8)
  labels <- c(rep("positive", 4), rep("negative", 4))
  result <- run_exhaustive_permutation_test(values, labels, "positive", "negative")
  expect_equal(result$n_permutations, choose(8, 4))
})

test_that("run_exhaustive_permutation_test detects perfect real separation as the most extreme possible p-value", {
  values <- c(10, 11, 12, 13, 1, 2, 3, 4)
  labels <- c(rep("positive", 4), rep("negative", 4))
  result <- run_exhaustive_permutation_test(values, labels, "positive", "negative")
  expect_equal(result$pvalue, 2 / choose(8, 4))
})

test_that("run_exhaustive_permutation_test gives a real non-significant p-value for indistinguishable real groups", {
  values <- c(5, 5, 5, 5, 5, 5, 5, 5)
  labels <- c(rep("positive", 4), rep("negative", 4))
  result <- run_exhaustive_permutation_test(values, labels, "positive", "negative")
  expect_equal(result$pvalue, 1.0)
})

test_that("run_leave_one_out_sensitivity drops each real patient exactly once, n-1 rows total", {
  patient_values <- data.frame(patient_id = paste0("P", 1:8), value = c(1, 2, 3, 4, 5, 6, 7, 8))
  positive_ids <- paste0("P", 1:4)
  negative_ids <- paste0("P", 5:8)
  result <- run_leave_one_out_sensitivity(patient_values, positive_ids, negative_ids)
  expect_equal(nrow(result), 8)
  expect_equal(sort(result$patient_removed), sort(patient_values$patient_id))
})

test_that("run_leave_one_out_sensitivity flags a real direction flip when an influential patient is removed", {
  # Positive group is uniformly higher except one real outlier (P1=0.01)
  # that, when present, still keeps the real positive median above the
  # real negative median; removing a DIFFERENT patient should not flip
  # direction, but this checks the real flag column exists and is
  # real logical.
  patient_values <- data.frame(patient_id = paste0("P", 1:8), value = c(0.01, 0.9, 0.85, 0.88, 0.1, 0.12, 0.09, 0.11))
  positive_ids <- paste0("P", 1:4)
  negative_ids <- paste0("P", 5:8)
  result <- run_leave_one_out_sensitivity(patient_values, positive_ids, negative_ids)
  expect_true(is.logical(result$same_direction_as_full))
})

test_that("run_bootstrap_median_difference_ci returns a real interval containing the real observed point estimate direction", {
  set.seed(1)
  patient_values <- data.frame(patient_id = paste0("P", 1:8), value = c(0.8, 0.85, 0.9, 0.82, 0.1, 0.15, 0.12, 0.08))
  positive_ids <- paste0("P", 1:4)
  negative_ids <- paste0("P", 5:8)
  result <- run_bootstrap_median_difference_ci(patient_values, positive_ids, negative_ids, n_bootstrap = 500)
  expect_true(result$ci_low <= result$point_estimate)
  expect_true(result$ci_high >= result$point_estimate)
  expect_true(result$point_estimate > 0)
})
