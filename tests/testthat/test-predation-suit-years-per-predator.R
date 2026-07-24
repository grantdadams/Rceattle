# Per-predator suitability reference-year windows (suit_styr / suit_endyr).
#
# These were scalar (a single global window) prior to v4.9.x; they are now a
# per-predator vector so each predator can average its suitability over a
# different set of years (e.g. hake 1980-2019, ATF 2013-2018, sablefish
# 2005-2008 in the California Current multispecies model). A scalar is recycled
# to every predator, which must reproduce the old behaviour exactly.

test_that("a recycled scalar suit window equals the legacy scalar behaviour", {
  dat <- make_msm_test_data()$data_list
  sy <- dat$styr
  ey <- dat$endyr
  np <- dat$nspp

  build <- function(ss, se) {
    fit_mod(
      data_list = dat, inits = NULL, file = NULL,
      estimateMode = 3, random_rec = FALSE, msmMode = 1,
      suitMode = 0, suit_styr = ss, suit_endyr = se,
      initMode = 2, niter = 3,
      fit_control = fit_control(verbose = 0)
    )
  }

  scalar_fit   <- build(sy, ey)
  recycled_fit <- build(rep(sy, np), rep(ey, np))

  # Scalar input is stored expanded to one value per predator ...
  expect_length(scalar_fit$data_list$suit_styr, np)
  expect_equal(scalar_fit$data_list$suit_styr, rep(sy, np))

  # ... and gives an identical objective to an explicit recycled vector.
  expect_equal(scalar_fit$obj$fn(), recycled_fit$obj$fn(), tolerance = 1e-12)
})

test_that("a distinct per-predator window changes the objective", {
  dat <- make_msm_test_data()$data_list
  sy <- dat$styr
  ey <- dat$endyr
  np <- dat$nspp
  mid <- as.integer((sy + ey) / 2)

  build <- function(ss, se) {
    fit_mod(
      data_list = dat, inits = NULL, file = NULL,
      estimateMode = 3, random_rec = FALSE, msmMode = 1,
      suitMode = 0, suit_styr = ss, suit_endyr = se,
      initMode = 2, niter = 3,
      fit_control = fit_control(verbose = 0)
    )
  }

  full_window <- build(rep(sy, np), rep(ey, np))
  # Narrow only the first predator's suitability-averaging window.
  narrow_pred1 <- build(c(mid, rep(sy, np - 1)), rep(ey, np))

  expect_equal(narrow_pred1$data_list$suit_styr, c(mid, rep(sy, np - 1)))
  expect_false(isTRUE(all.equal(full_window$obj$fn(),
                                narrow_pred1$obj$fn(),
                                tolerance = 1e-8)))
})

test_that("fit_mod rejects a per-predator window where styr > endyr", {
  dat <- make_msm_test_data()$data_list
  sy <- dat$styr
  ey <- dat$endyr
  np <- dat$nspp
  # Invert the first predator's window (styr > endyr); the rest are valid.
  expect_error(
    fit_mod(
      data_list = dat, inits = NULL, file = NULL,
      estimateMode = 3, random_rec = FALSE, msmMode = 1,
      suitMode = 0,
      suit_styr  = c(ey, rep(sy, np - 1)),
      suit_endyr = c(sy, rep(ey, np - 1)),
      initMode = 2, niter = 3,
      fit_control = fit_control(verbose = 0)
    ),
    "suit_styr must be <= suit_endyr"
  )
})

test_that("fit_mod rejects a window outside the hindcast [styr, endyr]", {
  dat <- make_msm_test_data()$data_list
  sy <- dat$styr
  ey <- dat$endyr
  np <- dat$nspp

  run <- function(ss, se) {
    fit_mod(
      data_list = dat, inits = NULL, file = NULL,
      estimateMode = 3, random_rec = FALSE, msmMode = 1,
      suitMode = 0, suit_styr = ss, suit_endyr = se,
      initMode = 2, niter = 3, fit_control = fit_control(verbose = 0))
  }

  # Start before styr would read the year arrays out of bounds (segfault in
  # MakeADFun) if it reached the C++; data_check must stop it first.
  expect_error(run(c(sy - 3, rep(sy, np - 1)), rep(ey, np)),
               "suit_styr must be >= styr")
  # End past endyr silently pads the averaging divisor with empty years.
  expect_error(run(rep(sy, np), c(ey + 5, rep(ey, np - 1))),
               "suit_endyr must be <= endyr")
})

test_that("fit_mod rejects a suit window vector of the wrong length", {
  dat <- make_msm_test_data()$data_list          # nspp == 2
  expect_error(
    fit_mod(
      data_list = dat, inits = NULL, file = NULL,
      estimateMode = 3, random_rec = FALSE, msmMode = 1,
      suitMode = 0,
      suit_styr  = rep(dat$styr, dat$nspp + 1),   # length 3 with nspp 2
      suit_endyr = rep(dat$endyr, dat$nspp),
      initMode = 2, niter = 3, fit_control = fit_control(verbose = 0)),
    "length-2 \\(nspp\\) vector; got length 3"
  )
})
