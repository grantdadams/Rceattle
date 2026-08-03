# =============================================================================
# M1_re = 6 (ar1_age_year) is a SEPARABLE 2D-AR1 on natural mortality: it must
# estimate BOTH an age correlation and a year correlation. A gate bug
# (build_map: `if(M1_re_model == 5)` inside the `%in% c(3, 6)` block, which can
# never be 5) plus an undefined index left mode 6's rho parameters unmapped, so
# the separable AR1 silently collapsed to IID (identical to mode 3). This pins
# that each AR1 M1_re mode frees the right number of rho hyperparameters:
#   mode 3 (iid_age_year)     -> 0 rho
#   mode 4 (ar1_age)          -> 1 rho
#   mode 5 (ar1_year)         -> 1 rho
#   mode 6 (ar1_age_year)     -> 2 rho   (the regression: was 0)
# =============================================================================

testthat::test_that("M1_re AR1 modes free the correct number of rho parameters", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  d <- make_test_data()
  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  npar <- function(mode) {
    f <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, M1Fun = Rceattle::build_M1(M1_model = 1, M1_re = mode),
      estimateMode = 3, msmMode = 0, random_rec = FALSE, fit_control = ctl)))
    testthat::expect_true(is.finite(f$obj$fn()))
    length(f$obj$par)
  }

  base <- npar(3)                                  # iid_age_year: 0 rho
  testthat::expect_equal(npar(4), base + 1)        # ar1_age:   +1 rho
  testthat::expect_equal(npar(5), base + 1)        # ar1_year:  +1 rho
  testthat::expect_equal(npar(6), base + 2)        # separable: +2 rho (was +0)
})
