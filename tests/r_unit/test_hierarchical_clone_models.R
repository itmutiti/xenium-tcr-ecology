# Unit tests for scripts/13_clone_ecology_confirmatory_models/02_fit_hierarchical_clone_models.R
# (`13_clone_ecology_confirmatory_models/02_fit_hierarchical_clone_models.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_hierarchical_clone_models.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/13_clone_ecology_confirmatory_models/02_fit_hierarchical_clone_models.R")

test_that("check_categorical_structure_confirmed is FALSE when all analyses are continuous", {
  df <- data.frame(analysis = c("primary", "sensitivity"), final_decision = c("continuous", "continuous"))
  expect_false(check_categorical_structure_confirmed(df))
})

test_that("check_categorical_structure_confirmed is TRUE only when ALL analyses are discrete", {
  df <- data.frame(analysis = c("primary", "sensitivity"), final_decision = c("discrete", "discrete"))
  expect_true(check_categorical_structure_confirmed(df))
})

test_that("check_categorical_structure_confirmed is FALSE on a real disagreement between analyses", {
  df <- data.frame(analysis = c("primary", "sensitivity"), final_decision = c("discrete", "continuous"))
  expect_false(check_categorical_structure_confirmed(df))
})
