# Unit tests for scripts/06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R
# (`06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_tme_substates.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R")

test_that("assign_substate_within_compartment picks the highest-scoring substate", {
  df <- data.frame(
    Myeloid__Macrophage_score = c(0.8, -0.2),
    Myeloid__Monocyte_score = c(0.1, 0.9)
  )
  result <- assign_substate_within_compartment(df, "Myeloid", c("Macrophage", "Monocyte"))
  expect_equal(result, c("Macrophage", "Monocyte"))
})

test_that("assign_substate_within_compartment handles 3 candidate substates", {
  df <- data.frame(
    Dendritic_cell__cDC_score = c(0.9, -0.5, 0.1),
    Dendritic_cell__pDC_score = c(0.1, 0.8, 0.1),
    Dendritic_cell__Mature_DC_score = c(-0.3, -0.1, 0.95)
  )
  result <- assign_substate_within_compartment(df, "Dendritic_cell", c("cDC", "pDC", "Mature_DC"))
  expect_equal(result, c("cDC", "pDC", "Mature_DC"))
})

test_that("ties break to the first-listed substate deterministically", {
  df <- data.frame(
    Endothelial__Blood_endothelial_score = c(0.5),
    Endothelial__Lymphatic_endothelial_score = c(0.5)
  )
  result <- assign_substate_within_compartment(df, "Endothelial", c("Blood_endothelial", "Lymphatic_endothelial"))
  expect_equal(result, "Blood_endothelial")
})

test_that("every compartment in COMPARTMENT_SUBSTATES has at least 2 candidate substates", {
  for (compartment in names(COMPARTMENT_SUBSTATES)) {
    expect_gte(length(COMPARTMENT_SUBSTATES[[compartment]]), 2)
  }
})

test_that("z-scoring recovers a genuine low-magnitude signal that raw-scale comparison would always miss (regression test)", {
  # Regression test for a real bug found on real data: a single-gene
  # score (e.g. Mature_DC from LAMP3 alone) is not on the same scale as a
  # multi-gene averaged score (e.g. cDC from 5 genes) -- scanpy's
  # score_genes does not shrink a 1-gene score toward zero the way
  # averaging does for a multi-gene score, so comparing raw scores directly
  # produced an implausible 75.7% "Mature_DC" call rate. Construct a case
  # where cDC's genuine signal (elevated, but only relative to its own
  # tight, low-variance baseline) is smaller in absolute magnitude than
  # Mature_DC's uniform offset -- raw-scale comparison must always pick
  # Mature_DC here (0% correct by construction), while z-scoring within
  # each column, which is exactly the point of the fix, must recover the
  # elevated cDC half using each column's own distribution.
  set.seed(1)
  n <- 200
  mature_dc_like <- rnorm(n, mean = 0.77, sd = 0.05)
  # cDC-like: a tight low baseline, and a "positive" half only modestly
  # elevated (mean 0.3, well below Mature_DC's 0.77) -- by construction,
  # raw-scale comparison can never pick cDC for either half.
  cdc_like <- c(rnorm(n / 2, mean = 0.05, sd = 0.02), rnorm(n / 2, mean = 0.3, sd = 0.02))

  # Raw (non-standardised) comparison: confirms the premise -- cDC can
  # never win here on a raw scale, exactly the bug's mechanism.
  raw_result <- ifelse(cdc_like > mature_dc_like, "cDC", "Mature_DC")
  expect_equal(mean(raw_result == "cDC"), 0)

  df <- data.frame(
    Dendritic_cell__cDC_score = cdc_like,
    Dendritic_cell__Mature_DC_score = mature_dc_like
  )
  result <- assign_substate_within_compartment(df, "Dendritic_cell", c("cDC", "Mature_DC"))

  # After z-scoring, the elevated cDC half (relative to its own
  # tight baseline) must be recoverable as cDC despite never winning on a
  # raw scale (0%, asserted above) -- a robust majority-recovery bar, not
  # 100%: standardisation equalises each column's mean/variance, not the
  # full distributional shape, so recovery is strong but imperfect by
  # construction, which is expected, not a sign the fix is incomplete.
  fraction_cdc_in_second_half <- mean(result[(n / 2 + 1):n] == "cDC")
  expect_gt(fraction_cdc_in_second_half, 0.75)
})
