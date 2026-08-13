# Unit tests for scripts/04_quality_control/06_define_qc_thresholds_hierarchically.R
# (`04_quality_control/06_define_qc_thresholds_hierarchically.R`). Run from the repository root:
#   Rscript -e 'testthat::test_file("tests/r_unit/test_qc_thresholds.R")'
# First R test file in the project -- see environment/conda/main.yml for the
# r-testthat justification comment (11+ scaffolded R scripts across later
# phases will follow this same convention).

suppressPackageStartupMessages(library(testthat))

# Sourcing (rather than running via Rscript) does not invoke main(): its
# call is guarded by `if (sys.nframe() == 0)`, which is only true for
# top-level Rscript execution, not for source(). Path is relative to this
# test file's own directory: testthat::test_file() sets the working
# directory to the test file's own directory (tests/r_unit/), not the
# directory Rscript was launched from -- confirmed empirically.
source("../../scripts/04_quality_control/06_define_qc_thresholds_hierarchically.R")

test_that("modified_z_scores flags a real outlier", {
  values <- c(10.0, 10.5, 9.8, 10.2, 9.9, 50.0)
  scores <- modified_z_scores(values)
  expect_gt(abs(scores[6]), 3.5)
  expect_true(all(abs(scores[1:5]) < 3.5))
})

test_that("modified_z_scores returns zero, not NaN/Inf, when MAD is zero", {
  values <- c(5.0, 5.0, 5.0, 5.0)
  scores <- modified_z_scores(values)
  expect_equal(scores, rep(0, 4))
})

test_that("modified_z_scores does not flag a tight cluster", {
  values <- c(10.0, 10.1, 9.9, 10.05, 9.95, 10.02)
  scores <- modified_z_scores(values)
  expect_true(all(abs(scores) < 3.5))
})

make_test_df <- function() {
  data.frame(
    transcript_counts = c(50, 3, 50, 50, 50),
    n_genes_detected = c(10, 10, 1, 10, 10),
    control_probe_ratio = c(0, 0, 0, 0.20, 0),
    control_codeword_ratio = c(0, 0, 0, 0, 0),
    z_counts = c(0, 0, 0, 0, -10)
  )
}

test_that("evaluate_profile flags low transcript_counts", {
  df <- make_test_df()
  excluded <- evaluate_profile(df, PROFILES$standard)
  expect_true(excluded[2])
})

test_that("evaluate_profile flags low n_genes_detected", {
  df <- make_test_df()
  excluded <- evaluate_profile(df, PROFILES$standard)
  expect_true(excluded[3])
})

test_that("evaluate_profile flags high control_probe_ratio", {
  df <- make_test_df()
  excluded <- evaluate_profile(df, PROFILES$standard)
  expect_true(excluded[4])
})

test_that("evaluate_profile flags a section-relative low-count outlier via z_counts", {
  df <- make_test_df()
  excluded <- evaluate_profile(df, PROFILES$standard)
  expect_true(excluded[5])
})

test_that("evaluate_profile does not flag a clean cell", {
  df <- make_test_df()
  excluded <- evaluate_profile(df, PROFILES$standard)
  expect_false(excluded[1])
})

test_that("strict profile excludes at least as many cells as standard, which excludes at least as many as lenient", {
  df <- make_test_df()
  n_lenient <- sum(evaluate_profile(df, PROFILES$lenient))
  n_standard <- sum(evaluate_profile(df, PROFILES$standard))
  n_strict <- sum(evaluate_profile(df, PROFILES$strict))
  expect_lte(n_lenient, n_standard)
  expect_lte(n_standard, n_strict)
})

test_that("write_yaml() output is valid YAML for a standard-compliant parser, even when a string value contains a backtick", {
  # Regression coverage for the failure class this project actually hit:
  # `config/qc_thresholds.yaml`'s `applied_at` field once held an unquoted,
  # backtick-wrapped markdown-style code span (a data value containing
  # display-formatting characters) that R's yaml::write_yaml() emitted
  # without quoting -- valid enough for R's own lenient yaml::read_yaml()
  # to read back, but rejected outright by a standard-compliant loader
  # (Python's PyYAML, the actual downstream consumer in
  # src/xenium_tcr_ecology/qc/apply_filters.py and friends). The real fix
  # was removing the backticks from the source value (they were
  # decorative, not meaningful data); this test additionally guards the
  # general pattern -- with the yaml package version installed when this
  # test was written (r-yaml 2.3.12), write_yaml() already auto-quotes
  # this specific case correctly, so this test is expected to pass today,
  # but exists to catch a recurrence if a future yaml package version (or
  # a differently-structured value) regresses the quoting behaviour again.
  tmp <- tempfile(fileext = ".yaml")
  on.exit(unlink(tmp))

  risky_config <- list(
    flags = list(
      example_flag = list(
        definition = "example",
        applied_at = "`04_quality_control/07_apply_qc_filters_with_audit_trail.py`"
      )
    )
  )
  write_yaml(risky_config, tmp)

  # R's own reader is lenient -- this alone would NOT have caught the bug.
  expect_no_error(read_yaml(tmp))

  # The real, stricter cross-language check: this is what actually broke.
  python3 <- Sys.which("python3")
  skip_if(python3 == "", "python3 not on PATH -- cannot run the cross-language YAML strictness check")
  result <- system2(
    python3,
    c("-c", shQuote("import sys, yaml; yaml.safe_load(open(sys.argv[1]))"), shQuote(tmp)),
    stdout = TRUE, stderr = TRUE
  )
  exit_code <- attr(result, "status")
  exit_code <- if (is.null(exit_code)) 0L else exit_code
  expect_equal(
    exit_code, 0L,
    info = paste(
      "config/qc_thresholds.yaml (or any file written the same way) must parse",
      "with a standard-compliant YAML loader, not just R's own yaml::read_yaml().",
      "python3 -c yaml.safe_load(...) reported:", paste(result, collapse = "\n")
    )
  )
})

test_that("the actual committed config/qc_thresholds.yaml parses with a standard-compliant YAML loader", {
  # Direct check against the real, tracked file -- not just a synthetic
  # reproduction of the bug pattern above. Skipped if the file hasn't been
  # generated yet in this environment (e.g. a from-scratch clean-room run
  # that hasn't reached phase04_06 yet); if present, it must be valid.
  qc_thresholds_path <- file.path(dirname(dirname(getwd())), "config", "qc_thresholds.yaml")
  candidate_paths <- c(qc_thresholds_path, file.path("..", "..", "config", "qc_thresholds.yaml"))
  existing <- candidate_paths[file.exists(candidate_paths)]
  skip_if(length(existing) == 0, "config/qc_thresholds.yaml does not exist in this environment yet")

  python3 <- Sys.which("python3")
  skip_if(python3 == "", "python3 not on PATH -- cannot run the cross-language YAML strictness check")
  result <- system2(
    python3,
    c("-c", shQuote("import sys, yaml; yaml.safe_load(open(sys.argv[1]))"), shQuote(existing[1])),
    stdout = TRUE, stderr = TRUE
  )
  exit_code <- attr(result, "status")
  exit_code <- if (is.null(exit_code)) 0L else exit_code
  expect_equal(
    exit_code, 0L,
    info = paste("config/qc_thresholds.yaml failed strict YAML parsing:", paste(result, collapse = "\n"))
  )
})
