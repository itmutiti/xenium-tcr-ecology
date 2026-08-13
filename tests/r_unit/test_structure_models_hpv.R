# Unit tests for scripts/15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R
# (`15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_structure_models_hpv.R")'

suppressPackageStartupMessages({
  library(testthat)
})

source("../../scripts/15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R")

test_that("aggregate_patient_category_metric averages equally across a patient's own sections", {
  data <- data.frame(
    section_id = c("S1", "S2", "S1", "S2"),
    category = c("A", "A", "B", "B"),
    value = c(0.2, 0.8, 0.4, 0.6)
  )
  section_to_patient <- c(S1 = "P1", S2 = "P1")
  result <- aggregate_patient_category_metric(data$section_id, data$category, data$value, section_to_patient)
  expect_equal(result$value[result$category == "A"], 0.5)
  expect_equal(result$value[result$category == "B"], 0.5)
})

test_that("aggregate_patient_clone_structure collapses section-then-patient (not raw clone-row mean)", {
  # Section S1 has 2 real clones with scores 0.0 and 1.0 (section mean 0.5).
  # Section S2 (same patient) has 1 real clone with score 0.9.
  # Real patient value should be equal-weight mean of the two SECTION means: (0.5+0.9)/2 = 0.7,
  # not the naive all-clone-row mean ((0+1+0.9)/3 = 0.633), which would let S1's extra clone dominate.
  data <- data.frame(
    section_id = c("S1", "S1", "S2"),
    ecological_structure_score = c(0.0, 1.0, 0.9)
  )
  section_to_patient <- c(S1 = "P1", S2 = "P1")
  result <- aggregate_patient_clone_structure(data$section_id, data$ecological_structure_score, section_to_patient)
  expect_equal(result$value[result$patient_id == "P1"], 0.7)
})

test_that("compare_hpv_groups detects a real engineered category-level difference", {
  patient_values <- data.frame(
    patient_id = c("Ppos1", "Ppos2", "Ppos3", "Ppos4", "Pneg1", "Pneg2", "Pneg3", "Pneg4"),
    category = "ecosystem_x",
    value = c(0.8, 0.85, 0.9, 0.82, 0.1, 0.15, 0.12, 0.08)
  )
  result <- compare_hpv_groups(patient_values, positive_ids = c("Ppos1", "Ppos2", "Ppos3", "Ppos4"), negative_ids = c("Pneg1", "Pneg2", "Pneg3", "Pneg4"))
  expect_true(result$median_positive > result$median_negative)
  expect_true(result$pvalue < 0.05)
})

test_that("compare_hpv_groups real p-value never falls below the real exact minimum for n=4 vs n=4", {
  patient_values <- data.frame(
    patient_id = c("Ppos1", "Ppos2", "Ppos3", "Ppos4", "Pneg1", "Pneg2", "Pneg3", "Pneg4"),
    category = "ecosystem_x",
    value = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)
  )
  result <- compare_hpv_groups(patient_values, positive_ids = c("Ppos1", "Ppos2", "Ppos3", "Ppos4"), negative_ids = c("Pneg1", "Pneg2", "Pneg3", "Pneg4"))
  expect_true(result$pvalue >= 2 / choose(8, 4) - 1e-9)
})
