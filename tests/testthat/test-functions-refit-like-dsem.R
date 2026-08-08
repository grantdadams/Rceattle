# .refit_like() carries the DSEM through every diagnostic refit.
#
# All eight refit entry points -- retrospective(), jitter(), self_test(),
# profile(), run_mse(), remove_F(), sample_rec(), reweight_comps() -- funnel
# through .refit_like(). Dropping the DSEM there would refit a DIFFERENT model
# than the one being diagnosed and report Mohn's rho / profiles / MSE output for
# it, with no error. These pin that it is forwarded, and forwarded as the
# SPECIFICATION rather than as built objects.

testthat::test_that(".refit_like() takes the DSEM spec from the data list by default", {
  # The default must read data_list$dsem_settings -- that is the whole mechanism
  # by which a refit inherits the DSEM without every caller passing it.
  fmls <- formals(Rceattle:::.refit_like)
  testthat::expect_true("dsem" %in% names(fmls))
  testthat::expect_identical(deparse(fmls$dsem), "data_list$dsem_settings")
})

testthat::test_that("the DSEM spec actually reaches fit_mod() through the refit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Behavioral, not textual: put a spec on the data list and check the refit
  # lands in fit_mod()'s DSEM block. Until the C++ is wired that block stops with
  # a known message, which is exactly the evidence we need -- reaching it means
  # the spec was forwarded and rebuilt against the refit's own data list.
  base <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))

  dl <- base$data_list
  dl$dsem_settings <- Rceattle::build_DSEM()   # default IID sem
  testthat::expect_error(
    suppressWarnings(suppressMessages(Rceattle:::.refit_like(
      data_list = dl, inits = base$estimated_params, estimateMode = 3))),
    "not yet wired into the compiled model")

  # ...and an explicit dsem = NULL overrides the data list, so a caller can opt
  # out. This is what makes the default a default rather than a hard-wire.
  testthat::expect_s3_class(
    suppressWarnings(suppressMessages(Rceattle:::.refit_like(
      data_list = dl, inits = base$estimated_params, estimateMode = 3,
      dsem = NULL))),
    "Rceattle")
})

testthat::test_that("a non-DSEM refit is unaffected", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # .refit_like() expects a data list already resolved by fit_mod() (scalar
  # srr_fun, HCR fields populated), which is how every real caller reaches it.
  base <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  dl <- base$data_list
  testthat::expect_null(dl$dsem_settings)

  fit <- suppressWarnings(suppressMessages(Rceattle:::.refit_like(
    data_list = dl, inits = base$estimated_params, estimateMode = 3)))
  testthat::expect_s3_class(fit, "Rceattle")
  testthat::expect_null(fit$dsem)
  # And it is the same model: mode 3 returns the real objective.
  testthat::expect_equal(fit$obj$fn(), base$obj$fn())
})

testthat::test_that("every refitting diagnostic refuses a DSEM with a directed message", {
  # Regression: these failed in assorted opaque ways instead of saying so --
  # retrospective() died zeroing inits$rec_dev (length 0 under a DSEM, since the
  # deviations live in the latent states), jitter() reported "0 of N returned"
  # with the cause buried, self_test() threw "argument is of length zero".
  fake <- structure(
    list(data_list = c(Rceattle::BS2017SS,
                       list(dsem_settings = Rceattle::build_DSEM(
                         sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1")))),
    class = "Rceattle")

  # A DSEM can live in three places; a guard that checks only dsem_settings is
  # silent on a fit built from pre-built objects, whose canonical slot is $dsem.
  built_only <- structure(
    list(data_list = Rceattle::BS2017SS,
         dsem = list(tmb_inputs = list(parameters = list(beta_z = 1)))),
    class = "Rceattle")
  testthat::expect_error(Rceattle::retrospective(built_only),
                         "does not yet support a DSEM")

  for (fn in c("retrospective", "jitter", "self_test", "remove_F")) {
    testthat::expect_error(do.call(fn, list(fake)), "does not yet support a DSEM",
                           info = fn)
  }
  testthat::expect_error(Rceattle::reweight_comps(fake), "does not yet support a DSEM")
  testthat::expect_error(stats::profile(fake, param = "sigmaR"),
                         "does not yet support a DSEM")
})

testthat::test_that("run_mse() and sample_rec() refuse a DSEM rather than mis-project it", {
  # Recruitment deviations are DERIVED from the latent states under a DSEM, so
  # sample_rec()'s draw into rec_dev would be overwritten on the next objective
  # evaluation -- the sampled recruitment silently unused. And whether x_tj is
  # projyr-dimensioned depends on build_DSEM(estimate_projection), which the
  # static .mse_proj_param_yrdim table cannot express. Both are deferred, so
  # both must fail loudly.
  fake <- structure(
    list(data_list = c(Rceattle::BS2017SS,
                       list(dsem_settings = Rceattle::build_DSEM(
                         sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1")))),
    class = "Rceattle")

  testthat::expect_error(Rceattle::run_mse(om = fake, em = fake, nsim = 1),
                         "does not yet support a DSEM")
  testthat::expect_error(Rceattle::sample_rec(fake),
                         "does not yet support a DSEM")
})
