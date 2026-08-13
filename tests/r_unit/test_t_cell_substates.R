# Unit tests for scripts/06_cell_type_annotation/04_resolve_t_cell_substates.R
# (`06_cell_type_annotation/04_resolve_t_cell_substates.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_t_cell_substates.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/06_cell_type_annotation/04_resolve_t_cell_substates.R")

make_row <- function(treg_score = -1, proliferation_score = -1, exhaustion_score = -1,
                      cytotoxicity_score = -1, cd4_expr = 0, cd8a_expr = 0) {
  data.frame(
    treg_score = treg_score, proliferation_score = proliferation_score,
    exhaustion_score = exhaustion_score, cytotoxicity_score = cytotoxicity_score,
    cd4_expr = cd4_expr, cd8a_expr = cd8a_expr
  )
}

test_that("a CD4-dominant, Treg-marker-positive cell is called Treg", {
  df <- make_row(treg_score = 0.5, cd4_expr = 2.0, cd8a_expr = 0.1)
  expect_equal(assign_t_cell_substate(df), "Treg")
})

test_that("a Treg-marker-positive but CD8-dominant cell is NOT called Treg", {
  # Tregs are conventionally CD4+ -- CD8-dominance should override a
  # positive Treg score.
  df <- make_row(treg_score = 0.5, cd4_expr = 0.1, cd8a_expr = 2.0)
  result <- assign_t_cell_substate(df)
  expect_false(result == "Treg")
  expect_equal(result, "CD8")
})

test_that("Cycling takes priority over Exhausted and Cytotoxic", {
  df <- make_row(proliferation_score = 0.5, exhaustion_score = 0.5, cytotoxicity_score = 0.5)
  expect_equal(assign_t_cell_substate(df), "Cycling")
})

test_that("Exhausted takes priority over Cytotoxic", {
  df <- make_row(exhaustion_score = 0.5, cytotoxicity_score = 0.5)
  expect_equal(assign_t_cell_substate(df), "Exhausted")
})

test_that("a cell with only a positive cytotoxicity score is called Cytotoxic", {
  df <- make_row(cytotoxicity_score = 0.5)
  expect_equal(assign_t_cell_substate(df), "Cytotoxic")
})

test_that("a cell with no positive program score falls back to CD4/CD8 by marker dominance", {
  df <- make_row(cd4_expr = 1.5, cd8a_expr = 0.2)
  expect_equal(assign_t_cell_substate(df), "CD4")

  df2 <- make_row(cd4_expr = 0.2, cd8a_expr = 1.5)
  expect_equal(assign_t_cell_substate(df2), "CD8")
})

test_that("a cell with no positive scores and no CD4/CD8 expression is Ambiguous", {
  df <- make_row()
  expect_equal(assign_t_cell_substate(df), "Ambiguous")
})

test_that("assign_t_cell_substate is vectorised across multiple rows", {
  df <- rbind(
    make_row(treg_score = 0.5, cd4_expr = 2.0, cd8a_expr = 0.1),
    make_row(cytotoxicity_score = 0.5),
    make_row()
  )
  result <- assign_t_cell_substate(df)
  expect_equal(result, c("Treg", "Cytotoxic", "Ambiguous"))
})
