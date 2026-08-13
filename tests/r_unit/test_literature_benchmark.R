# Unit tests for scripts/14_spatial_interactions_and_barriers/06_benchmark_against_published_barrier_studies.R
# (`14_spatial_interactions_and_barriers/06_benchmark_against_published_barrier_studies.R`, `q3_literature_benchmark` -- supporting, not independently
# multiplicity-corrected). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_literature_benchmark.R")'

suppressPackageStartupMessages({
  library(testthat)
})

source("../../scripts/14_spatial_interactions_and_barriers/06_benchmark_against_published_barrier_studies.R")

test_that("compute_raw_correlation recovers a known real engineered correlation", {
  set.seed(1)
  x <- rnorm(200)
  y <- -0.5 * x + rnorm(200, sd = 0.5)
  result <- compute_raw_correlation(x, y)
  expect_true(result$r < -0.3)
  expect_true(result$r > -0.9)
  expect_true(result$pvalue < 0.05)
})

test_that("compute_raw_correlation recovers near-zero r for real independent noise", {
  set.seed(2)
  x <- rnorm(200)
  y <- rnorm(200)
  result <- compute_raw_correlation(x, y)
  expect_true(abs(result$r) < 0.2)
})

test_that("classify_concordance flags real same-sign correlations as direction-concordant", {
  result <- classify_concordance(project_r = -0.08, published_r = -0.48)
  expect_equal(result$direction_concordant, TRUE)
  expect_equal(result$magnitude_ratio, 0.08 / 0.48)
})

test_that("classify_concordance flags real opposite-sign correlations as direction-discordant", {
  result <- classify_concordance(project_r = 0.08, published_r = -0.48)
  expect_equal(result$direction_concordant, FALSE)
})

test_that("classify_concordance real magnitude_ratio is NA-safe when published_r is zero", {
  result <- classify_concordance(project_r = 0.1, published_r = 0)
  expect_true(is.na(result$magnitude_ratio))
})
