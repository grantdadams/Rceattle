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
  # covr instruments the TMB model at -O0, which historically moved the GOA fits
  # to a different point of their flat selectivity ridge. Polishing to a
  # stationary point ought to make the optima build-independent, but that is not
  # verified here, so the skip is kept.
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Pinned objective functions (this branch), at Newton-POLISHED optima.
  # Without newtonsteps each fit stops on nlminb's objective-RELATIVE tolerance,
  # so it halts wherever that tolerance happens to bite rather than at a
  # stationary point -- the GOA fits in particular sat ~1e-3 in gradient, far
  # above the package's own 1e-4 convergence threshold. That made the reference
  # sensitive to changes that cannot alter the model at all: adding a constant to
  # the objective (which leaves every gradient untouched) once moved goa_ss by
  # 52.9 units. Polishing pins a true optimum instead, which is reproducible
  # under any such perturbation. getsd = FALSE is numerically inert here.
  ref <- c(ss     = 10241.0304272585,
           ms     = 10267.2478324443,
           goa_ss = 12868.0052289274,
           goa_ms = 12932.7931701136)
  tol <- 1e-6

  fc <- function(...) Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                            newtonsteps = 3, ...)

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

  fits <- list(ss = ss, ms = ms, goa_ss = goa_ss, goa_ms = goa_ms)
  got  <- vapply(fits, function(f) f$opt$objective, numeric(1))

  for (m in names(ref))
    testthat::expect_equal(got[[m]], ref[[m]], tolerance = tol, info = m)

  # Each constant must pin a CONVERGED optimum, not merely a reproducible number.
  # Without this the pinned values can drift back to a non-stationary point and
  # still pass, which is how the reference silently became fragile before.
  for (m in names(ref))
    testthat::expect_lt(max(abs(fits[[m]]$obj$gr(fits[[m]]$opt$par))), 1e-4)
})
