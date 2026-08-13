# Unit tests for scripts/11_clone_spatial_descriptors/06_discover_provisional_structure.R
# (`11_clone_spatial_descriptors/06_discover_provisional_structure.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_provisional_structure.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/11_clone_spatial_descriptors/06_discover_provisional_structure.R")

test_that("compute_continuous_factor_scores returns one real score per observation", {
  set.seed(1)
  loadings <- matrix(rnorm(6), ncol = 1)
  factor_scores_true <- rnorm(100)
  data <- outer(factor_scores_true, as.vector(loadings)) + matrix(rnorm(100 * 6, sd = 0.2), ncol = 6)
  result <- compute_continuous_factor_scores(data, n_factors = 1)
  expect_equal(length(result), 100)
  expect_true(all(is.finite(result)))
})

test_that("compute_continuous_factor_scores recovers the real underlying latent axis", {
  set.seed(2)
  loadings <- matrix(c(1, 1, 1, 1, 1, 1), ncol = 1)
  factor_scores_true <- rnorm(200)
  data <- outer(factor_scores_true, as.vector(loadings)) + matrix(rnorm(200 * 6, sd = 0.1), ncol = 6)
  result <- compute_continuous_factor_scores(data, n_factors = 1)
  # Recovered scores should correlate strongly (up to sign) with the real
  # true latent axis used to generate the data.
  correlation <- abs(cor(result, factor_scores_true))
  expect_true(correlation > 0.9)
})

test_that("extract_factor_loadings returns one row per feature, ordered by |loading|", {
  set.seed(3)
  loadings_true <- matrix(c(2, 0.1, 1.5, 0.2, 1, 0.05), ncol = 1)
  factor_scores_true <- rnorm(150)
  data <- outer(factor_scores_true, as.vector(loadings_true)) + matrix(rnorm(150 * 6, sd = 0.1), ncol = 6)
  colnames(data) <- paste0("feature_", 1:6)
  result <- extract_factor_loadings(data, n_factors = 1)
  expect_equal(nrow(result), 6)
  expect_true(all(diff(abs(result$loading)) <= 1e-8))
  # The largest real loading (feature_1, true loading 2) should rank first.
  expect_equal(result$feature[1], "feature_1")
})
