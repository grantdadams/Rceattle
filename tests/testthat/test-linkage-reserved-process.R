# Regression tests for the reserved-but-unwired "sel" process.
#
# "sel" is reserved in LINKAGE_PROCESSES (and as RCEATTLE_PROC_SEL in
# src/TMB/linkage.hpp) but has no accumulator yet. A "sel" row must be
# rejected up front rather than estimated without affecting the model, or
# dying with an opaque "unknown param" deeper in the pipeline.
#
# ("q" was in the same state until catchability linkages were wired; it is
# now implemented and has its own tests in test-linkage-catchability.R.)

.env_df <- function() data.frame(Year = 1:6, temp = seq(-1, 1, length.out = 6))

testthat::test_that("sel is reserved but not yet implemented", {
  testthat::expect_true("sel" %in% Rceattle:::LINKAGE_PROCESSES)
  testthat::expect_false("sel" %in% Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED)
  # q has been wired; sel has not.
  testthat::expect_true("q" %in% Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED)
  testthat::expect_setequal(Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED,
                            c("recruitment", "M", "growth", "q"))
})


testthat::test_that("materialize_linkage() rejects the sel process", {
  spec <- Rceattle::linkage_spec(~temp, param = "M1")
  testthat::expect_error(
    Rceattle:::materialize_linkage(spec, "sel", .env_df(),
                                   list(species = 1L)),
    "reserved but not yet implemented"
  )
})


testthat::test_that("materialize_linkage() still works for wired processes", {
  spec <- Rceattle::linkage_spec(~temp, param = "M1")
  tbl <- Rceattle:::materialize_linkage(spec, "M", .env_df(),
                                        list(species = 1L))
  testthat::expect_s3_class(tbl, "Rceattle_linkage_table")
  # intercept + slope
  testthat::expect_equal(nrow(tbl), 2L)
  testthat::expect_true(all(tbl$process == "M"))
})


testthat::test_that("validate_linkage_table() backstops hand-built sel rows", {
  testthat::expect_error(
    Rceattle:::linkage_row(process = "sel", param = "coff", X_col = 1L),
    "reserved but not yet implemented"
  )
})
