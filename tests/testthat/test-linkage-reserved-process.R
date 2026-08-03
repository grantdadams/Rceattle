# Every reserved linkage process is now wired to a C++ accumulator, so the
# reserved-but-unimplemented guard should currently reject nothing. These tests
# pin that state and the two guards that remain meaningful:
#   * a genuinely unknown process still errors, and
#   * an unimplemented *link* function still errors.
# (recruitment/M/growth/q have their own tests; sel is in
# test-linkage-selectivity.R.)

.env_df <- function() data.frame(Year = 1:6, temp = seq(-1, 1, length.out = 6))

testthat::test_that("all reserved processes are now implemented", {
  # LINKAGE_PROCESSES is the reserved set; every one now has an accumulator,
  # so IMPLEMENTED equals the full reserved set.
  testthat::expect_setequal(Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED,
                            Rceattle:::LINKAGE_PROCESSES)
  testthat::expect_true(all(c("recruitment", "M", "growth", "q", "sel") %in%
                              Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED))
})


testthat::test_that("an unknown process is still rejected", {
  # A process outside the reserved set must never encode.
  testthat::expect_error(
    Rceattle:::linkage_row(process = "not_a_process", param = "x", X_col = 1L))
})


testthat::test_that("materialize_linkage() works for a wired process", {
  spec <- Rceattle::linkage_spec(~temp, param = "M1")
  tbl <- Rceattle:::materialize_linkage(spec, "M", .env_df(),
                                        list(species = 1L))
  testthat::expect_s3_class(tbl, "Rceattle_linkage_table")
  testthat::expect_equal(nrow(tbl), 2L)   # intercept + slope
  testthat::expect_true(all(tbl$process == "M"))
})


testthat::test_that("the logit link is still reserved (no accumulator)", {
  # Every accumulator gates on log (1) / identity (0); logit (2) has no
  # consumer, so a logit-link row must be rejected rather than silently
  # contribute zero while its base parameter is mapped out.
  testthat::expect_error(
    Rceattle:::linkage_row(process = "M", param = "M1", X_col = 1L,
                           link = "logit"),
    "reserved|not yet implemented|logit")
})
