# Unit tests for scripts/16_external_validation_and_generalisation/03_validate_cell_state_signatures.R
# (`16_external_validation_and_generalisation/03_validate_cell_state_signatures.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_signature_validation.R")'

suppressPackageStartupMessages({
  library(testthat)
  library(Matrix)
})

source("../../scripts/16_external_validation_and_generalisation/03_validate_cell_state_signatures.R")

test_that("identify_t_cells_from_marker_counts flags cells with any real positive marker count", {
  # 3 real marker genes (rows) x 4 real cells (columns).
  marker_counts <- Matrix(c(
    0, 0, 0, 0,
    1, 0, 0, 0,
    0, 0, 2, 0
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  result <- identify_t_cells_from_marker_counts(marker_counts)
  expect_equal(result, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("normalize_counts_cpm_log scales by real per-cell library size then log1p", {
  counts <- c(10, 20, 0)
  total_counts <- c(100, 200, 50)
  median_total <- 100
  result <- normalize_counts_cpm_log(counts, total_counts, median_total)
  # cell 1: 10/100*100 = 10 -> log1p(10); cell 2: 20/200*100 = 10 -> log1p(10) (same normalised value despite different raw counts/depth)
  expect_equal(result[1], log1p(10))
  expect_equal(result[2], log1p(10))
  expect_equal(result[3], log1p(0))
})

test_that("compute_mean_pairwise_correlation recovers a known real value on perfectly correlated genes", {
  set.seed(1)
  base <- rnorm(50)
  expr_matrix <- cbind(base, base + rnorm(50, sd = 0.001), base + rnorm(50, sd = 0.001))
  result <- compute_mean_pairwise_correlation(expr_matrix)
  expect_true(result > 0.99)
})

test_that("compute_mean_pairwise_correlation gives a near-zero value for real independent genes", {
  set.seed(2)
  expr_matrix <- matrix(rnorm(50 * 5), ncol = 5)
  result <- compute_mean_pairwise_correlation(expr_matrix)
  expect_true(abs(result) < 0.3)
})

test_that("compute_module_coherence_pvalue detects a real engineered coherent module against a real incoherent background", {
  set.seed(3)
  n_cells <- 100
  base <- rnorm(n_cells)
  # Real observed module: 4 real, highly correlated genes.
  observed_expr <- cbind(base, base + rnorm(n_cells, sd = 0.05), base + rnorm(n_cells, sd = 0.05), base + rnorm(n_cells, sd = 0.05))
  observed_corr <- compute_mean_pairwise_correlation(observed_expr)
  # Real background: 200 real, mutually independent genes.
  background_expr <- matrix(rnorm(n_cells * 200), ncol = 200)
  result <- compute_module_coherence_pvalue(observed_corr, background_expr, n_genes = 4, n_permutations = 200)
  expect_true(result$pvalue < 0.05)
})
