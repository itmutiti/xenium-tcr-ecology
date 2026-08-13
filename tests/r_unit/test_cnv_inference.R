# Unit tests for scripts/07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R
# (`07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_cnv_inference.R")'

suppressPackageStartupMessages(library(testthat))

source("../../scripts/07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R")

test_that("genes_by_chromosome splits genes by chromosome, excludes chrY", {
  gene_coords <- data.frame(
    gene = c("A", "B", "C", "D"),
    chromosome = c("1", "1", "2", "Y"),
    stringsAsFactors = FALSE
  )
  result <- genes_by_chromosome(gene_coords)
  expect_equal(sort(result[["1"]]), c("A", "B"))
  expect_equal(result[["2"]], "C")
  expect_false("Y" %in% names(result))
})

test_that("compute_relative_expression subtracts each cell's own patient's baseline", {
  epi_expr <- data.frame(
    patient_id = c("P1", "P1", "P2"),
    GENE_A = c(5, 7, 10),
    GENE_B = c(1, 2, 3),
    row.names = c("cellA", "cellB", "cellC")
  )
  reference_baseline <- data.frame(GENE_A = c(4, 8), GENE_B = c(0, 1), row.names = c("P1", "P2"))
  result <- compute_relative_expression(epi_expr, reference_baseline, c("GENE_A", "GENE_B"))
  expect_equal(unname(result["cellA", "GENE_A"]), 1)
  expect_equal(unname(result["cellB", "GENE_A"]), 3)
  expect_equal(unname(result["cellC", "GENE_A"]), 2)
  expect_equal(unname(result["cellC", "GENE_B"]), 2)
})

test_that("compute_relative_expression errors on an unmatched patient_id", {
  epi_expr <- data.frame(patient_id = "P99", GENE_A = 1, row.names = "cellX")
  reference_baseline <- data.frame(GENE_A = 4, row.names = "P1")
  expect_error(compute_relative_expression(epi_expr, reference_baseline, "GENE_A"), "no matching")
})

test_that("compute_chromosome_cnv_scores averages within chromosome", {
  relative_expression <- matrix(
    c(1, 3, 5, 2, 4, 6),
    nrow = 2, dimnames = list(c("c1", "c2"), c("A", "B", "C"))
  )
  chrom_gene_list <- list("1" = c("A", "B"), "2" = "C")
  result <- compute_chromosome_cnv_scores(relative_expression, chrom_gene_list)
  expect_equal(unname(result["c1", "1"]), mean(c(1, 5)))
  expect_equal(unname(result["c1", "2"]), 4)
})

test_that("compute_cnv_burden is the variance across chromosomes, higher for more divergent cells", {
  chromosome_scores <- matrix(
    c(0, 0, 0, -5, 0, 5),
    nrow = 2, byrow = TRUE, dimnames = list(c("flat_cell", "divergent_cell"), c("1", "2", "3"))
  )
  result <- compute_cnv_burden(chromosome_scores)
  expect_equal(unname(result["flat_cell"]), 0)
  expect_gt(unname(result["divergent_cell"]), 0)
})
