# =============================================================================
# Integer-coded model switches accept readable string aliases (the CIE review's
# "rename options to improve interpretability, e.g. constant -> not estimated").
# switch_check() resolves the strings to the integer codes the model consumes,
# so a string spec must produce a BIT-IDENTICAL fit to the integer spec. New
# aliases use a consistent convention: the not-estimated value is "Fixed",
# estimated is "Estimated".
# =============================================================================

testthat::test_that("estDynamics resolves readable strings to integer codes", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  for (pair in list(c("Estimated", "0"), c("Fixed", "1"),
                    c("FixedScaled", "2"), c("FixedScaledByAge", "3"))) {
    d$estDynamics <- pair[1]
    out <- suppressMessages(Rceattle::switch_check(d))
    testthat::expect_equal(out$estDynamics, as.integer(pair[2]), info = pair[1])
  }
})

testthat::test_that("an unknown estDynamics string errors clearly", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  d$estDynamics <- "Frozen"
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "Invalid 'estDynamics'")
})

testthat::test_that("estDynamics string spec fits bit-identically to the integer spec", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  d_int <- make_test_data(seed = 3)
  d_int$estDynamics <- 0L
  d_str <- d_int
  d_str$estDynamics <- "Estimated"

  f_int <- suppressMessages(Rceattle::fit_mod(
    data_list = d_int, estimateMode = 3, msmMode = 0, random_rec = FALSE, fit_control = ctl))
  f_str <- suppressMessages(Rceattle::fit_mod(
    data_list = d_str, estimateMode = 3, msmMode = 0, random_rec = FALSE, fit_control = ctl))
  testthat::expect_equal(sum(f_str$quantities$jnll_comp),
                         sum(f_int$quantities$jnll_comp), tolerance = 0)
})
