# Unit tests for scripts/13_clone_ecology_confirmatory_models/03_test_clone_size_as_confounder.R
# (`13_clone_ecology_confirmatory_models/03_test_clone_size_as_confounder.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_clone_size_confounder.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/13_clone_ecology_confirmatory_models/03_test_clone_size_as_confounder.R")

test_that("compute_proportion_shift is zero when nothing changes", {
  props <- list(patient_proportion = 0.3, identity_proportion = 0.2, context_proportion = 0.5)
  result <- compute_proportion_shift(props, props)
  expect_equal(result$patient_shift, 0)
  expect_equal(result$identity_shift, 0)
  expect_equal(result$context_shift, 0)
})

test_that("compute_proportion_shift reflects a real, engineered shift", {
  unadjusted <- list(patient_proportion = 0.3, identity_proportion = 0.2, context_proportion = 0.5)
  adjusted <- list(patient_proportion = 0.1, identity_proportion = 0.2, context_proportion = 0.7)
  result <- compute_proportion_shift(unadjusted, adjusted)
  expect_equal(result$patient_shift, -0.2)
  expect_equal(result$identity_shift, 0)
  expect_equal(result$context_shift, 0.2)
})

test_that("fit_size_adjusted_model includes a real log_n_cells fixed effect", {
  set.seed(1)
  # Real replication: each of 20 clones appears twice (matching this
  # project's own real technical-replicate structure) -- a clone-level
  # random effect is not identifiable with only one row per clone.
  n_clones <- 20
  data <- data.frame(
    patient_id = rep(paste0("P", rep(1:5, each = 4)), each = 2),
    clone_id = rep(paste0("C", 1:n_clones), each = 2),
    n_cells = sample(1:50, n_clones * 2, replace = TRUE),
    ecological_structure_score = rnorm(n_clones * 2)
  )
  model <- fit_size_adjusted_model(data)
  expect_true("log_n_cells" %in% names(lme4::fixef(model)))
})
