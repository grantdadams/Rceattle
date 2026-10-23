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


testthat::test_that("a start with no positive innovation variance is refused", {
  testthat::skip_if_not_installed("dsem")

  # Both take a square root of the innovation variance, sigma^2 (1 - rho^2) at
  # lag 1: 0 at rho = 1, where the build aborts.
  d <- .cv_data()
  sem_at <- function(rho, name = "rhoR", extra = "") paste0(
    "recdevs1 -> recdevs1, 1, ", name, ", ", rho, "\n", extra,
    "recdevs1 <-> recdevs1, 0, sigmaR1, 1")
  build <- function(sem, cv) suppressWarnings(build_dsem_objects(
    build_DSEM(sem = sem, family = "fixed", constant_variance = cv),
    data_list = d))
  msg <- "no positive innovation variance holds the marginal variance"

  testthat::expect_error(build(sem_at(1), "marginal"), msg)
  testthat::expect_error(build(sem_at(-1), "diagonal"), msg)
  testthat::expect_error(build(sem_at(1.2), "marginal"), msg)
  testthat::expect_error(build(sem_at(1, "NA"), "marginal"), msg)
  # Each lag below 1, but together they do not damp.
  lag2 <- function(r) paste0("recdevs1 -> recdevs1, 2, rho2, ", r, "\n")
  testthat::expect_error(build(sem_at(0.6, extra = lag2(0.6)), "diagonal"), msg)
  # Stationary, or the conditional reading, which has no marginal to hold.
  testthat::expect_length(build(sem_at(0.9), "marginal")$sem_full$path, 2L)
  testthat::expect_length(build(sem_at(0.3, extra = lag2(0.2)), "marginal")$sem_full$path, 3L)
  testthat::expect_length(build(sem_at(1), "conditional")$sem_full$path, 2L)
})


testthat::test_that("a warm start at |rho| >= 1 is refused before the model is built", {
  testthat::skip_if_not_installed("dsem")

  # The sem starts at 0.9, so build_dsem_objects() passes; the inits are what
  # the tape would be built with.
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  d <- e$make_test_data(nyrs = 12L, nprojyrs = 3L, nages = 8L, minage = 1L,
                        growth = "vonBertalanffy", seed = 1)
  d$random_rec <- TRUE
  run <- function(inits = NULL) suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, inits = inits,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    dsem = build_DSEM(sem = paste0("recdevs1 -> recdevs1, 1, rhoR, 0.9\n",
                                   "recdevs1 <-> recdevs1, 0, sigmaR1, 1"),
                      constant_variance = "marginal"),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))
  fit <- run()
  ip  <- fit$estimated_params
  k   <- fit$dsem$sem_full$parameter[fit$dsem$sem_full$name %in% "rhoR"]
  testthat::expect_equal(ip$dsem_beta_z[k], 0.9)
  # The fit's own inits build; the same inits at rho = 1 do not.
  testthat::expect_true(is.list(run(ip)$obj))
  ip$dsem_beta_z[k] <- 1
  testthat::expect_error(run(ip),
                         "no positive innovation variance holds the marginal variance")
})


# Provenance: "marginal" failed on these shapes with bad_alloc in getParameterOrder.
testthat::test_that("marginal and diagonal hold each marginal variance on a two-lag sem", {
  testthat::skip_if_not_installed("dsem")

  # Each latent keeps its `<->` variance in every year under both settings.
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  d <- e$make_test_data(nyrs = 12L, nprojyrs = 3L, nages = 8L, minage = 1L,
                        growth = "vonBertalanffy", seed = 1)
  d$random_rec <- TRUE
  fit_with <- function(sem, cv) suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    dsem = build_DSEM(sem = sem, constant_variance = cv),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))
  two <- paste0(
    "recdevs1 <-> recdevs1, 0, sR1, 1\nrecdevs1 -> recdevs1, 1, rho, 0.6\n",
    "mdevs1 <-> mdevs1, 0, sM, 0.2\nmdevs1 -> mdevs1, 1, rM, 0.6\n")
  for (cv in c("marginal", "diagonal")) {
    fit <- fit_with(two, cv)
    dt  <- fit$dsem$tmb_inputs$data
    mv  <- fit$obj$report(fit$obj$env$par)$dsem_margvar_tj
    testthat::expect_equal(as.numeric(mv[, dt$rec_dev_col[1] + 1L]),
                           rep(1, nrow(mv)), tolerance = 1e-8, info = cv)
    testthat::expect_equal(as.numeric(mv[, dt$mdev_col + 1L]),
                           rep(0.04, nrow(mv)), tolerance = 1e-8, info = cv)
  }
  lag2 <- fit_with("recdevs1 <-> recdevs1, 0, sR1, 1\nrecdevs1 -> recdevs1, 2, rho2, 0.3\n",
                   "marginal")
  testthat::expect_true(is.finite(lag2$obj$fn(lag2$obj$par)))

  # Past the boundary a double evaluation reports NaN, not a stand-in model.
  one <- fit_with("recdevs1 <-> recdevs1, 0, sR1, 1\nrecdevs1 -> recdevs1, 1, rho, 0.6\n",
                  "marginal")
  p <- one$obj$env$par
  k <- which(names(p) == "dsem_beta_z")
  testthat::expect_length(k, 1L)
  testthat::expect_true(is.finite(one$obj$report(p)$jnll))
  p[k] <- 1
  r <- one$obj$report(p)
  testthat::expect_true(is.na(r$jnll))
  testthat::expect_true(all(is.na(r$R[1, seq_len(one$dsem$tmb_inputs$data$nyrs_dsem)])))

  # A covariance between two variables: "marginal" would scale it wrongly.
  cov <- paste0("recdevs1 <-> recdevs1, 0, sR1, 1\nmdevs1 <-> mdevs1, 0, sM, 0.2\n",
                "recdevs1 <-> mdevs1, 0, cRM, 0.3\n")
  testthat::expect_error(fit_with(cov, "marginal"), "mis-scales a covariance")
  dg <- fit_with(cov, "diagonal")
  mv <- dg$obj$report(dg$obj$env$par)$dsem_margvar_tj
  testthat::expect_equal(mv[1, dg$dsem$tmb_inputs$data$rec_dev_col[1] + 1L], 1.09,
                         tolerance = 1e-8)
})
