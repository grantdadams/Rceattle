# =============================================================================
# Golden-reference regression: the four reference fits (Bering Sea + Gulf of
# Alaska, single- and multi-species) must reproduce their pinned objective
# functions. This is the committed backstop for "numeric changes need golden
# equivalence" -- any edit that moves a reference fit fails here.
#
# The multispecies fits are warm-started from their single-species MLEs (the
# predation likelihood is non-convex, so the start point selects the local
# optimum); the GOA multispecies fit uses fixed M. See the /golden-check
# command for the recipe these constants were generated from.
# =============================================================================

testthat::test_that("the four reference fits reproduce their pinned objectives", {
  testthat::skip_on_cran()
  # covr instruments the TMB model at -O0; the GOA fits then land on a different
  # point of their flat selectivity ridge, so the pinned -O2 optima don't hold.
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Pinned objective functions (this branch). getsd = FALSE is numerically inert
  # for the objective and much faster.
  ref <- c(ss     = 10241.0304275402,
           ms     = 10267.2478327352,
           goa_ss = 12807.4375258732,
           goa_ms = 12866.2957391829)
  tol <- 1e-6

  fc <- function(...) Rceattle::fit_control(getsd = FALSE, verbose = 0, ...)

  ss <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, file = NULL,
    inits = NULL, estimateMode = 0, random_rec = FALSE, msmMode = 0,
    fit_control = fc(phase = TRUE))
  ms <- Rceattle::fit_mod(data_list = Rceattle::BS2017MS,
    inits = ss$estimated_params, file = NULL, estimateMode = 0, niter = 5,
    random_rec = FALSE, msmMode = 1, suitMode = 0, fit_control = fc())
  goa_ss <- Rceattle::fit_mod(data_list = Rceattle::GOA2018SS, file = NULL,
    inits = NULL, estimateMode = 0, random_rec = FALSE, msmMode = 0,
    fit_control = fc(phase = TRUE))
  goa_ms <- Rceattle::fit_mod(data_list = Rceattle::GOA2018SS,
    inits = goa_ss$estimated_params, file = NULL, estimateMode = 0, niter = 3,
    random_rec = FALSE, msmMode = 1, suitMode = 0, fit_control = fc(phase = TRUE))

  got <- c(ss = ss$opt$objective, ms = ms$opt$objective,
           goa_ss = goa_ss$opt$objective, goa_ms = goa_ms$opt$objective)

  for (m in names(ref))
    testthat::expect_equal(got[[m]], ref[[m]], tolerance = tol, info = m)
})
