# Unit tests for scripts/15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R
# (`15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_composition_models.R")'

suppressPackageStartupMessages({
  library(testthat)
})

source("../../scripts/15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R")

test_that("assign_section_id_from_cell_id matches the real unique section_id prefix", {
  cell_ids <- c("P01_run1_aaadggoi-1", "P09_run2_zzzzzz-1")
  valid_sections <- c("P01_run1", "P09_run1", "P09_run2")
  result <- assign_section_id_from_cell_id(cell_ids, valid_sections)
  expect_equal(result, c("P01_run1", "P09_run2"))
})

test_that("assign_section_id_from_cell_id errors on an unmatched cell_id", {
  cell_ids <- c("P99_run1_zzzzzz-1")
  valid_sections <- c("P01_run1", "P09_run1")
  expect_error(assign_section_id_from_cell_id(cell_ids, valid_sections))
})

test_that("compute_section_lineage_fractions sums to 1.0 per section", {
  lineage <- c("T_cell", "T_cell", "Fibroblast", "T_cell", "Fibroblast")
  section_id <- c("S1", "S1", "S1", "S2", "S2")
  result <- compute_section_lineage_fractions(lineage, section_id, all_lineages = c("T_cell", "Fibroblast"))
  s1 <- result[result$section_id == "S1", ]
  expect_equal(s1$fraction[s1$lineage == "T_cell"], 2 / 3)
  expect_equal(s1$fraction[s1$lineage == "Fibroblast"], 1 / 3)
  totals <- tapply(result$fraction, result$section_id, sum)
  expect_equal(as.numeric(totals), c(1.0, 1.0))
})

test_that("aggregate_patient_lineage_fractions averages equally across a patient's sections (not cell-count-weighted)", {
  section_fractions <- data.frame(
    section_id = c("S1", "S2"),
    lineage = c("T_cell", "T_cell"),
    fraction = c(0.2, 0.8)
  )
  section_to_patient <- c(S1 = "P1", S2 = "P1")
  result <- aggregate_patient_lineage_fractions(section_fractions, section_to_patient)
  expect_equal(result$fraction[result$patient_id == "P1" & result$lineage == "T_cell"], 0.5)
})

test_that("compare_hpv_groups_by_lineage detects a real engineered group difference", {
  set.seed(1)
  patient_fractions <- data.frame(
    patient_id = rep(c("Ppos1", "Ppos2", "Ppos3", "Ppos4", "Pneg1", "Pneg2", "Pneg3", "Pneg4"), each = 1),
    lineage = "T_cell",
    fraction = c(0.8, 0.85, 0.9, 0.82, 0.1, 0.15, 0.12, 0.08)
  )
  result <- compare_hpv_groups_by_lineage(patient_fractions, positive_ids = c("Ppos1", "Ppos2", "Ppos3", "Ppos4"), negative_ids = c("Pneg1", "Pneg2", "Pneg3", "Pneg4"))
  row <- result[result$lineage == "T_cell", ]
  expect_true(row$median_positive > row$median_negative)
  expect_true(row$pvalue < 0.05)
})

test_that("compare_hpv_groups_by_lineage real p-values are never below the real exact minimum for n=4 vs n=4", {
  patient_fractions <- data.frame(
    patient_id = c("Ppos1", "Ppos2", "Ppos3", "Ppos4", "Pneg1", "Pneg2", "Pneg3", "Pneg4"),
    lineage = "T_cell",
    fraction = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)
  )
  result <- compare_hpv_groups_by_lineage(patient_fractions, positive_ids = c("Ppos1", "Ppos2", "Ppos3", "Ppos4"), negative_ids = c("Pneg1", "Pneg2", "Pneg3", "Pneg4"))
  expect_true(result$pvalue[result$lineage == "T_cell"] >= 2 / choose(8, 4) - 1e-9)
})
