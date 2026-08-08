# fit_mod(dsem = ) plumbing.
#
# The C++ side has not landed, so a DSEM fit cannot run yet. What these pin is
# the half that is real now: that `dsem = NULL` is a complete no-op (the gate
# that lets golden-check police this change), that a DSEM warm start does not
# poison a non-DSEM fit, that supplying a `dsem` builds its inputs and then
# fails with a directed message rather than an opaque TMB one, and that the
# specification round-trips through the run configuration.

sem_for <- function(nspp) {
  paste(sprintf("recdevs%d <-> recdevs%d, 0, sigmaR%d, 1",
                seq_len(nspp), seq_len(nspp), seq_len(nspp)), collapse = "\n")
}

testthat::test_that("dsem defaults to NULL and leaves the objective untouched", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  testthat::expect_true("dsem" %in% names(formals(Rceattle::fit_mod)))
  testthat::expect_null(formals(Rceattle::fit_mod)$dsem)

  fc <- Rceattle::fit_control(getsd = FALSE, verbose = 0)
  base <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fc)))
  explicit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0, dsem = NULL,
    fit_control = fc)))

  # estimateMode 3 returns the real objective, so this compares the actual model.
  testthat::expect_identical(base$obj$fn(), explicit$obj$fn())
})

testthat::test_that("a DSEM warm start is scrubbed when the fit has no DSEM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Regression: DSEM parameters are absent from build_params(), so fit_mod()'s
  # inits shape guard passes them through as unflagged extras and MakeADFun then
  # rejects parameters the template never declares. Reached whenever inits come
  # from a DSEM fit (e.g. one produced on the dev-DSEM branch).
  fc <- Rceattle::fit_control(getsd = FALSE, verbose = 0)
  base <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fc)))

  poisoned <- base$estimated_params
  poisoned$beta_z    <- c(0.5, 0.8)
  poisoned$lnsigma_z <- 0
  poisoned$mu_j      <- c(0, 0)
  poisoned$delta0_j  <- c(0, 0)
  poisoned$x_tj      <- matrix(0, 5, 2)

  refit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = poisoned, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fc)))

  testthat::expect_false(any(c("beta_z", "lnsigma_z", "mu_j", "delta0_j", "x_tj")
                             %in% names(refit$initial_params)))
  testthat::expect_identical(refit$obj$fn(), base$obj$fn())
})

testthat::test_that("supplying a dsem builds its inputs, then reports the C++ gap", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  nyr <- d$endyr - d$styr + 1
  d$env_data <- data.frame(Year = d$styr:d$endyr,
                           BT = as.numeric(scale(seq_len(nyr))))

  # Reaching this error means build_dsem_objects() succeeded on a v5.0 data
  # list -- which is the point of the check, since the DSEM builder read the
  # pre-v5.0 `sigma_rec_prior` name and used to fail before getting here.
  testthat::expect_error(
    suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3,
      random_rec = TRUE, msmMode = 0,
      dsem = Rceattle::build_DSEM(sem = sem_for(d$nspp)),
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    "not yet wired into the compiled model"
  )
})

testthat::test_that("the DSEM spec round-trips through a run configuration", {
  spec <- Rceattle::build_DSEM(sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                               sigmaR_prior_sd = 0.5)
  cfg  <- Rceattle::model_config(msmMode = 0, dsem = spec)
  testthat::expect_false(is.null(cfg$dsem))

  f <- withr::local_tempfile(fileext = ".yaml")
  Rceattle::save_config(Rceattle::run_config(cfg), f)
  back <- Rceattle::load_config(f)

  testthat::expect_false(is.null(back$model_config$dsem))
  testthat::expect_identical(back$model_config$dsem$sem, spec$sem)
  testthat::expect_equal(back$model_config$dsem$sigmaR_prior_sd, 0.5)

  # An omitted dsem must not reach the written YAML, so a stored config can
  # never silently turn DSEM on. (model_config() keeps the slot as an explicit
  # NULL; what matters is that serialization drops it.)
  testthat::expect_null(Rceattle::model_config(msmMode = 0)$dsem)
  f2 <- withr::local_tempfile(fileext = ".yaml")
  Rceattle::save_config(Rceattle::run_config(Rceattle::model_config(msmMode = 0)), f2)
  testthat::expect_false(any(grepl("^dsem:", readLines(f2))))
  testthat::expect_null(Rceattle::load_config(f2)$model_config$dsem)
})

testthat::test_that("a DSEM and a recruitment linkage are rejected together", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # Regression: this compared linkage_table$process (a process NAME) against
  # LINKAGE_PROCESS_CODES[["recruitment"]] (integer 0), which coerces to "0" and
  # never matches -- so the check silently never fired. Once the C++ lands that
  # would let a DSEM and a recruitment linkage both structure the recruitment
  # deviations in the same fit.
  d <- Rceattle::BS2017SS
  nyr <- d$endyr - d$styr + 1
  d$env_data <- data.frame(Year = d$styr:d$endyr,
                           BT = as.numeric(scale(seq_len(nyr))))

  testthat::expect_error(
    suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3,
      random_rec = TRUE, msmMode = 0,
      recFun = Rceattle::build_srr(
        srr_fun  = 0,
        linkages = list(R0 = Rceattle::linkage_spec(~ BT, by = ~ species))),
      dsem = Rceattle::build_DSEM(sem = sem_for(d$nspp)),
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    "cannot be combined"
  )
})

testthat::test_that("dsem must be a build_DSEM() specification", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Regression: the type check used inherits(dsem, "Rceattle_dsem"), a class no
  # constructor ever assigns, so it was always FALSE -- any object was accepted
  # and a bare list with a $sem element had its family / estimate_projection
  # silently defaulted.
  testthat::expect_error(
    suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
      estimateMode = 3, random_rec = FALSE, msmMode = 0,
      dsem = "not a spec",
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    "must be a specification from build_DSEM"
  )
})

testthat::test_that("a diagnostic refit inherits the DSEM from the data list", {
  # .refit_like() rebuilds every spec from the data list; the DSEM must be among
  # them, or retrospective()/jitter()/profile() etc. would refit without it and
  # report diagnostics for a different model. Covered in depth in
  # test-functions-refit-like-dsem.R; this pins the wiring from this side.
  fmls <- formals(Rceattle:::.refit_like)
  testthat::expect_true("dsem" %in% names(fmls))
  testthat::expect_identical(deparse(fmls$dsem), "data_list$dsem_settings")
})

testthat::test_that("build_dsem_objects accepts the v5.0 sigma_rec spelling", {
  testthat::skip_if_not_installed("dsem")

  # Regression: v5.0 renamed the control element sigma_rec_prior -> sigma_rec
  # (schema alias, upgraded on entry), but the DSEM builder still read the old
  # name, so it saw NULL and rep_len() failed on every v5.0 data list.
  base <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L,
               spnames = "a", random_rec = TRUE, proj_mean_rec = TRUE,
               env_data = data.frame(Year = 1990:2000, BT = 0))
  sem <- "recdevs1 <-> recdevs1, 0, sigmaR1, 1"

  d_new <- base; d_new$sigma_rec <- 1
  d_old <- base; d_old$sigma_rec_prior <- 1
  d_neither <- base

  o_new <- build_dsem_objects(dsem_settings = build_DSEM(sem = sem), data_list = d_new)
  o_old <- build_dsem_objects(dsem_settings = build_DSEM(sem = sem), data_list = d_old)
  testthat::expect_equal(o_new$tmb_inputs$data$rec_sd_prior,
                         o_old$tmb_inputs$data$rec_sd_prior)

  testthat::expect_error(
    build_dsem_objects(dsem_settings = build_DSEM(sem = sem), data_list = d_neither),
    "neither `sigma_rec` nor `sigma_rec_prior`")
})
