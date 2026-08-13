# Unit tests for scripts/10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R
# (`10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_neighbourhood_archetypes.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R")

test_that("compute_consensus_matrix gives consensus 1 for well-separated clusters", {
  # Two tight, well-separated blobs -- every bootstrap resample should put
  # each point back with its own blob, giving near-perfect within-blob
  # consensus and near-zero between-blob consensus.
  set.seed(1)
  blob_a <- matrix(rnorm(20 * 2, mean = 0, sd = 0.05), ncol = 2)
  blob_b <- matrix(rnorm(20 * 2, mean = 10, sd = 0.05), ncol = 2)
  data <- rbind(blob_a, blob_b)

  consensus <- compute_consensus_matrix(data, k = 2, n_bootstrap = 20, subsample_frac = 0.8, seed = 1)

  within_a <- consensus[1:20, 1:20][upper.tri(consensus[1:20, 1:20])]
  between <- consensus[1:20, 21:40]
  expect_true(mean(within_a, na.rm = TRUE) > 0.95)
  expect_true(mean(between, na.rm = TRUE) < 0.05)
})

test_that("compute_consensus_matrix has NA diagonal and is symmetric", {
  set.seed(2)
  data <- matrix(rnorm(30 * 2), ncol = 2)
  consensus <- compute_consensus_matrix(data, k = 2, n_bootstrap = 10, seed = 2)
  expect_true(all(is.na(diag(consensus))))
  offdiag <- consensus[upper.tri(consensus)]
  expect_true(all(is.na(offdiag) | (offdiag >= 0 & offdiag <= 1)))
  expect_equal(consensus[1, 2], consensus[2, 1])
})

test_that("compute_pac_score is near zero for a perfectly bimodal consensus matrix", {
  # A consensus matrix with only 0s and 1s off-diagonal has no entries in
  # the ambiguous [0.1, 0.9] zone at all.
  consensus <- matrix(c(NA, 1, 0, 1, NA, 0, 0, 0, NA), nrow = 3)
  expect_equal(compute_pac_score(consensus), 0)
})

test_that("compute_pac_score is 1 when every entry is fully ambiguous", {
  consensus <- matrix(c(NA, 0.5, 0.5, 0.5, NA, 0.5, 0.5, 0.5, NA), nrow = 3)
  expect_equal(compute_pac_score(consensus), 1)
})

test_that("compute_pac_score ignores NA (never-co-sampled) entries", {
  consensus <- matrix(c(NA, NA, 1, NA, NA, 0, 1, 0, NA), nrow = 3)
  # Off-diagonal upper-tri entries: [1,2]=NA, [1,3]=1, [2,3]=0 -- only two
  # real observations, both outside the ambiguous zone.
  expect_equal(compute_pac_score(consensus), 0)
})

test_that("select_stable_k picks the K with minimum PAC", {
  pac_by_k <- c("4" = 0.3, "6" = 0.1, "8" = 0.25, "10" = 0.4)
  expect_equal(select_stable_k(pac_by_k), 6)
})

test_that("label_archetypes_by_dominant_lineage picks the max-fraction column", {
  centroids <- rbind(
    c(B_cell = 0.7, T_cell = 0.3),
    c(B_cell = 0.2, T_cell = 0.8)
  )
  labels <- label_archetypes_by_dominant_lineage(centroids, c("B_cell", "T_cell"))
  expect_equal(labels, c("B_cell", "T_cell"))
})
