# Regression tests for the reserved-but-unwired "q" and "sel" processes.
#
# Both are reserved in LINKAGE_PROCESSES (and as RCEATTLE_PROC_Q /
# RCEATTLE_PROC_SEL in src/TMB/linkage.hpp) but neither has an
# accumulator. They used to fail asymmetrically and both failed badly:
#   * LINKAGE_PARAM_CODES$q is populated, so a "q" row encoded cleanly
#     and was then estimated and prior-penalized while contributing
#     nothing to the model.
#   * LINKAGE_PARAM_CODES$sel is empty, so a "sel" row died inside
#     encode_linkage_for_tmb() with an opaque "unknown param".
# Both must now be rejected up front with a message that says why.

.env_df <- function() data.frame(Year = 1:6, temp = seq(-1, 1, length.out = 6))

testthat::test_that("q and sel are reserved but not implemented", {
  testthat::expect_true(all(c("q", "sel") %in% Rceattle:::LINKAGE_PROCESSES))
  testthat::expect_false(any(c("q", "sel") %in%
                               Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED))
  testthat::expect_setequal(Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED,
                            c("recruitment", "M", "growth"))
})


testthat::test_that("materialize_linkage() rejects reserved processes", {
  spec <- Rceattle::linkage_spec(~temp, param = "M1")
  for (proc in c("q", "sel")) {
    testthat::expect_error(
      Rceattle:::materialize_linkage(spec, proc, .env_df(),
                                     list(species = 1L)),
      "reserved but not yet implemented",
      info = proc
    )
  }
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


testthat::test_that("validate_linkage_table() backstops hand-built q/sel rows", {
  testthat::expect_error(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L),
    "reserved but not yet implemented"
  )
})
