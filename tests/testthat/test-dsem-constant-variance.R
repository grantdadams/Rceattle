# build_DSEM(constant_variance = ) -- what a `<->` term means once the sem has
# lagged paths.
#
# dsem inherits RAM notation from the `sem` package, where a `<->` self-loop is
# the DISTURBANCE variance, conditional on the incoming paths. So the marginal
# variance grows along the causal graph: sigma^2 / (1 - rho^2) for a first-order
# self-path. "marginal" / "diagonal" rescale Gamma and (I - Rho) to hold the
# marginal variance constant instead -- the time-series convention that
# ar1(1 | Year) already uses here, and WHAM's Ecov AR(1) and glmmTMB's ar1() use
# elsewhere.
#
# Rceattle keeps dsem's default. These tests pin (a) that the default has not
# moved, so no existing fit changes, (b) that the argument actually reaches the
# template's options vector rather than being accepted and dropped, (c) dsem's
# own claim that the three are equivalent without lags, and (d) that the setting
# survives a run-config round trip -- a build_DSEM() argument that save_config()
# does not carry round-trips to nothing, which is how index_cov was lost.

.cv_data <- function() {
  list(styr = 1990L, endyr = 2010L, projyr = 2015L, nspp = 1L, spnames = "a",
       sigma_rec = 1, random_rec = TRUE, proj_mean_rec = TRUE,
       env_data = data.frame(Year = 1990:2015, BT = as.numeric(1:26)))
}

testthat::test_that("the default is dsem's, so no existing fit moves", {
  testthat::expect_equal(build_DSEM()$constant_variance, "conditional")
  testthat::expect_equal(eval(formals(build_DSEM)$constant_variance)[1],
                         "conditional")
  # Same set dsem offers, in dsem's order.
  testthat::expect_equal(eval(formals(build_DSEM)$constant_variance),
                         c("conditional", "marginal", "diagonal"))
  testthat::expect_error(build_DSEM(constant_variance = "nope"))
})

testthat::test_that("the setting reaches the template's options vector", {
  testthat::skip_if_not_installed("dsem")
  d <- .cv_data()
  sem <- "BT -> BT, 1, AR_BT, 0
BT -> recdevs1, 1, BT_to_R, 0
recdevs1 <-> recdevs1, 0, sigmaR1, 1"
  # options(1) -> 0: constant conditional variance; 1: constant marginal; 2: diagonal
  for (i in seq_along(c("conditional", "marginal", "diagonal"))) {
    cv <- c("conditional", "marginal", "diagonal")[i]
    b <- suppressWarnings(build_dsem_objects(
      build_DSEM(sem = sem, family = "fixed", constant_variance = cv),
      data_list = d))
    testthat::expect_equal(b$tmb_inputs$data$options[2], i - 1L, info = cv)
  }

  # A spec built before this argument existed -- or assembled by hand -- must
  # fall back to dsem's default rather than reaching match.arg() with NULL.
  legacy <- list(sem = sem, family = "fixed", sigmaR_prior_sd = NA,
                 estimate_projection = FALSE)
  b <- suppressWarnings(build_dsem_objects(legacy, data_list = d))
  testthat::expect_equal(b$tmb_inputs$data$options[2], 0L)
})

testthat::test_that("the setting survives a run-config round trip", {
  spec <- build_DSEM(sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                     family = "fixed", constant_variance = "marginal")
  rc <- .rce_run_config(mc = model_config(dsem = spec))
  back <- .rce_run_config_from_list(.rce_run_config_to_list(rc))
  testthat::expect_equal(back$model_config$dsem$constant_variance, "marginal")
})

testthat::test_that("the three settings agree when the sem has no lags", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not_installed("TMB")

  # dsem's own documented claim: "equivalent when the model includes no lags
  # (only simultaneous effects) and no covariances".
  d <- Rceattle::GOApollock
  d$projyr <- 2020
  iid <- "recdevs1 <-> recdevs1, 0, sigmaR1, 0.7"
  obj <- vapply(c("conditional", "marginal", "diagonal"), function(cv) {
    m <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = "DebugBuild",
      random_rec = TRUE, msmMode = 0, initMode = 2,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0),
      dsem = Rceattle::build_DSEM(sem = iid, family = "fixed",
                                  constant_variance = cv))))
    as.numeric(m$obj$fn())
  }, numeric(1))
  testthat::expect_equal(obj[["marginal"]], obj[["conditional"]])
  testthat::expect_equal(obj[["diagonal"]],  obj[["conditional"]])
})

testthat::test_that("a lagged sem gives the marginal variance the convention implies", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not_installed("TMB")

  # rho is FIXED at 0.6 by the sem's start value and DebugBuild does not
  # optimize, so the closed form is exact rather than approximate.
  d <- Rceattle::GOApollock
  d$projyr <- 2020
  rho <- 0.6
  lagged <- paste0("recdevs1 -> recdevs1, 1, rhoR, ", rho, "\n",
                   "recdevs1 <-> recdevs1, 0, sigmaR1, 0.7")
  ratio <- vapply(c("conditional", "marginal"), function(cv) {
    m <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = "DebugBuild",
      random_rec = TRUE, msmMode = 0, initMode = 2,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0, bias_adjust_proc = TRUE),
      dsem = Rceattle::build_DSEM(sem = lagged, family = "fixed",
                                  constant_variance = cv))))
    col <- m$dsem$tmb_inputs$data$rec_dev_col[1] + 1L
    mv  <- as.numeric(m$quantities$dsem_margvar_tj[, col])
    mv[length(mv)] / as.numeric(m$quantities$R_sd)[1]^2
  }, numeric(1))

  # conditional: the <-> term is the innovation SD, so the marginal variance is
  # inflated by exactly 1/(1 - rho^2).
  testthat::expect_equal(ratio[["conditional"]], 1 / (1 - rho^2),
                         tolerance = 1e-6)
  # marginal: the <-> term IS the marginal SD, so the ratio is 1.
  testthat::expect_equal(ratio[["marginal"]], 1, tolerance = 1e-6)
})
