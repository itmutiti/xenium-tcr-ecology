# Unit tests for scripts/04_quality_control/08_assess_replicate_concordance.R (`04_quality_control/08_assess_replicate_concordance.R`).
# Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_replicate_concordance.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/04_quality_control/08_assess_replicate_concordance.R")

test_that("as_bool handles pandas-style 'True'/'False' character strings", {
  expect_equal(as_bool(c("True", "False", "True")), c(TRUE, FALSE, TRUE))
})

test_that("as_bool handles the factor representation arrow::read_parquet returns", {
  f <- factor(c("True", "False", "True"))
  expect_equal(as_bool(f), c(TRUE, FALSE, TRUE))
})

test_that("as_bool passes through an already-logical vector unchanged", {
  expect_equal(as_bool(c(TRUE, FALSE)), c(TRUE, FALSE))
})

test_that("jaccard_index is 1.0 for identical sets", {
  expect_equal(jaccard_index(c("a", "b", "c"), c("a", "b", "c")), 1.0)
})

test_that("jaccard_index is 0.0 for disjoint sets", {
  expect_equal(jaccard_index(c("a", "b"), c("c", "d")), 0.0)
})

test_that("jaccard_index computes the expected ratio for a partial overlap", {
  # intersection {b, c} = 2, union {a, b, c, d} = 4
  expect_equal(jaccard_index(c("a", "b", "c"), c("b", "c", "d")), 0.5)
})

test_that("jaccard_index handles two empty sets without dividing by zero", {
  expect_true(is.na(jaccard_index(character(0), character(0))))
})

test_that("modified_z_scores flags a real outlier", {
  values <- c(0.95, 0.96, 0.94, 0.97, 0.93, 0.40)
  scores <- modified_z_scores(values)
  expect_gt(abs(scores[6]), 3.5)
  expect_true(all(abs(scores[1:5]) < 3.5))
})

test_that("modified_z_scores returns zero, not NaN, when MAD is zero", {
  expect_equal(modified_z_scores(rep(0.9, 4)), rep(0, 4))
})

test_that("hull_area computes the area of a unit square exactly", {
  x <- c(0, 1, 1, 0)
  y <- c(0, 0, 1, 1)
  expect_equal(hull_area(x, y), 1.0)
})

test_that("hull_area computes the area of a right triangle exactly", {
  x <- c(0, 4, 0)
  y <- c(0, 0, 3)
  expect_equal(hull_area(x, y), 6.0)
})

test_that("moran_i returns a high positive value for perfectly spatially clustered values", {
  # Two well-separated clusters, each internally constant -- strong positive
  # spatial autocorrelation is the correct, expected result.
  coords <- rbind(
    cbind(runif(20, 0, 1), runif(20, 0, 1)),
    cbind(runif(20, 100, 101), runif(20, 100, 101))
  )
  values <- c(rep(0, 20), rep(100, 20))
  result <- moran_i(coords, values, k = 5)
  expect_gt(result, 0.5)
})

test_that("moran_i returns a value near zero for spatially random values", {
  set.seed(1)
  coords <- cbind(runif(200, 0, 100), runif(200, 0, 100))
  values <- rnorm(200)
  result <- moran_i(coords, values, k = 6)
  expect_lt(abs(result), 0.3)
})
