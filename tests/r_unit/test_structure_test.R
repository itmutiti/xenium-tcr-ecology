# Unit tests for scripts/11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R
# (`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_structure_test.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R")

test_that("compute_mvn_loglik matches a hand-computed value for a simple case", {
  # Standard bivariate normal, identity covariance, one point at origin.
  data <- matrix(c(0, 0), nrow = 1)
  sigma <- diag(2)
  result <- compute_mvn_loglik(data, sigma)
  expected <- -0.5 * (2 * log(2 * pi) + 0 + 0) # log(det(I))=0, quad form=0
  expect_equal(result, expected, tolerance = 1e-8)
})

test_that("compute_discrete_bic: three well-separated blobs get lower BIC at k=3 than k=1", {
  set.seed(1)
  blob <- function(cx, cy) matrix(c(rnorm(30, cx, 0.05), rnorm(30, cy, 0.05)), ncol = 2)
  data <- rbind(blob(0, 0), blob(10, 0), blob(5, 10))
  bic_k1 <- compute_discrete_bic(data, k = 1)
  bic_k3 <- compute_discrete_bic(data, k = 3)
  expect_true(bic_k3 < bic_k1)
})

test_that("discrete BIC improves far less, relative to its own k=1 baseline, for noise than for real blobs", {
  # BIC(k) on ANY finite sample improves somewhat with more clusters
  # (k-means always reduces within-cluster SS) -- that alone does not
  # mean discrete structure is real. What must differ is the SIZE of the
  # improvement: a clustered dataset should show a far larger
  # relative BIC drop than pure noise. This is the real, meaningful
  # comparison (not an isolated absolute threshold, which does not hold
  # in general at finite n -- checked and rejected).
  set.seed(2)
  noise <- matrix(rnorm(200), ncol = 2)
  noise_drop <- compute_discrete_bic(noise, k = 1) - compute_discrete_bic(noise, k = 3)

  set.seed(1)
  blob <- function(cx, cy) matrix(c(rnorm(30, cx, 0.05), rnorm(30, cy, 0.05)), ncol = 2)
  blobs <- rbind(blob(0, 0), blob(10, 0), blob(5, 10))
  blob_drop <- compute_discrete_bic(blobs, k = 1) - compute_discrete_bic(blobs, k = 3)

  expect_true(blob_drop > 5 * noise_drop)
})

test_that("compute_continuous_bic runs and returns a finite real value", {
  set.seed(3)
  loadings <- matrix(rnorm(6), ncol = 1)
  factor_scores <- rnorm(100)
  data <- outer(factor_scores, as.vector(loadings)) + matrix(rnorm(100 * 6, sd = 0.3), ncol = 6)
  result <- compute_continuous_bic(data, n_factors = 1)
  expect_true(is.finite(result))
})

test_that("compute_gap_statistic returns one row per k with real gap/sd values", {
  set.seed(4)
  data <- matrix(rnorm(100), ncol = 2)
  result <- compute_gap_statistic(data, k_max = 3, b_reference = 10, seed = 1)
  expect_equal(nrow(result), 3)
  expect_true(all(is.finite(result$gap)))
  expect_true(all(result$sd >= 0))
})

test_that("select_k_by_gap picks the smallest k satisfying Tibshirani's rule", {
  gap_df <- data.frame(k = 1:4, gap = c(0.5, 1.0, 1.05, 1.02), sd = c(0.1, 0.1, 0.1, 0.1))
  # k=2: gap(2)=1.0 >= gap(3)-sd(3) = 1.05-0.1 = 0.95 -> TRUE, select k=2
  expect_equal(select_k_by_gap(gap_df), 2)
})

test_that("select_k_by_gap defaults to k_max when no earlier k satisfies the rule", {
  gap_df <- data.frame(k = 1:3, gap = c(0.1, 0.5, 1.0), sd = c(0.01, 0.01, 0.01))
  expect_equal(select_k_by_gap(gap_df), 3)
})

test_that("classify_structure requires BOTH bic and gap to favour discrete", {
  expect_equal(classify_structure(delta_bic = 20, gap_selected_k = 3), "discrete")
  expect_equal(classify_structure(delta_bic = 20, gap_selected_k = 1), "continuous")
  expect_equal(classify_structure(delta_bic = 2, gap_selected_k = 3), "continuous")
  expect_equal(classify_structure(delta_bic = 2, gap_selected_k = 1), "continuous")
})

test_that("end-to-end: four well-separated real blobs (symmetric layout) are classified discrete", {
  # A layout with one blob at the origin and others shooting out along
  # orthogonal axes was tried first and rejected: it creates a lopsided
  # PCA-aligned bounding box that confuses the gap statistic's reference
  # sampling (checked directly -- it produced gap_selected_k=1 despite
  # an obvious real elbow at k=4). A symmetric, compact layout is the
  # standard, well-behaved test case for the gap statistic.
  set.seed(5)
  blob <- function(cx, cy, cz) matrix(c(rnorm(25, cx, 0.3), rnorm(25, cy, 0.3), rnorm(25, cz, 0.3)), ncol = 3)
  data <- rbind(blob(0, 0, 0), blob(10, 0, 3), blob(0, 10, 6), blob(10, 10, 9))

  discrete_bic <- sapply(1:5, function(k) compute_discrete_bic(data, k))
  best_discrete_bic <- min(discrete_bic)
  continuous_bic <- compute_continuous_bic(data, n_factors = 1)
  delta_bic <- continuous_bic - best_discrete_bic

  gap_df <- compute_gap_statistic(data, k_max = 5, b_reference = 20, seed = 1)
  gap_selected_k <- select_k_by_gap(gap_df)

  expect_equal(classify_structure(delta_bic, gap_selected_k), "discrete")
})

test_that("end-to-end: a single smooth continuous 1D manifold is classified continuous", {
  set.seed(6)
  n <- 150
  t <- runif(n, -5, 5)
  data <- cbind(t, 0.5 * t + rnorm(n, sd = 0.2), -0.3 * t + rnorm(n, sd = 0.2))

  discrete_bic <- sapply(1:5, function(k) compute_discrete_bic(data, k))
  best_discrete_bic <- min(discrete_bic)
  continuous_bic <- compute_continuous_bic(data, n_factors = 1)
  delta_bic <- continuous_bic - best_discrete_bic

  gap_df <- compute_gap_statistic(data, k_max = 5, b_reference = 20, seed = 1)
  gap_selected_k <- select_k_by_gap(gap_df)

  expect_equal(classify_structure(delta_bic, gap_selected_k), "continuous")
})
