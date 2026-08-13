# Unit tests for scripts/13_clone_ecology_confirmatory_models/04_run_leave_one_patient_out_stability.R
# (`13_clone_ecology_confirmatory_models/04_run_leave_one_patient_out_stability.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_lopo_stability.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/13_clone_ecology_confirmatory_models/04_run_leave_one_patient_out_stability.R")

test_that("all_components_non_trivial is TRUE when every component clears the threshold", {
  props <- list(patient_proportion = 0.3, identity_proportion = 0.2, context_proportion = 0.5)
  expect_true(all_components_non_trivial(props, threshold = 0.05))
})

test_that("all_components_non_trivial is FALSE when one real component falls below threshold", {
  props <- list(patient_proportion = 0.01, identity_proportion = 0.2, context_proportion = 0.79)
  expect_false(all_components_non_trivial(props, threshold = 0.05))
})

test_that("all_components_non_trivial is exactly TRUE at the boundary (>=)", {
  props <- list(patient_proportion = 0.05, identity_proportion = 0.3, context_proportion = 0.65)
  expect_true(all_components_non_trivial(props, threshold = 0.05))
})
