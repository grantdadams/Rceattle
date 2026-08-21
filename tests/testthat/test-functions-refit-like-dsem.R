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
  dl$env_data <- data.frame(Year = dl$styr:dl$endyr, BT = 0)
  dl$dsem_settings <- Rceattle::build_DSEM()   # default IID sem
  refit <- suppressWarnings(suppressMessages(Rceattle:::.refit_like(
    data_list = dl, inits = base$estimated_params, estimateMode = 3)))
  # The refit rebuilt the DSEM from the spec against its own data list. This is
  # the case that used to fail with "Map and parameter objects are not the same
  # size": `inits` from a DSEM fit already carry the dsem_* blocks, so appending
  # the template added a second copy of each.
  testthat::expect_s3_class(refit, "Rceattle")
  testthat::expect_false(is.null(refit$dsem))
  testthat::expect_gt(sum(refit$quantities$jnll_comp[22L, ]), 0)

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

testthat::test_that("the refitting diagnostics that cannot handle a DSEM refuse with a directed message", {
  # Regression: these failed in assorted opaque ways instead of saying so --
  # jitter() reported "0 of N returned" with the cause buried, self_test() threw
  # "argument is of length zero". They still refuse; retrospective() no longer
  # does, because its peel now marginalizes the latent states.
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
  # Pointed at a function that still refuses. remove_F() used to serve here and
  # no longer does, at which point this ran remove_F()'s body against a fake
  # object and failed on `max(NULL)` instead of testing the guard.
  testthat::expect_error(Rceattle::process_residuals(built_only),
                         "does not support a DSEM")
  # ...and process = "catchability" is NOT refused: index_q_dev is untouched by
  # a DSEM, so blocking it would cost a valid diagnostic and explain it with a
  # message about a process the caller did not ask for.
  testthat::expect_error(
    Rceattle::process_residuals(built_only, process = "catchability"),
    "no TMB object")

  # retrospective(), jitter() and self_test() are deliberately NOT in this list
  # any more. retrospective() peels by MARGINALIZING the peeled-year latent
  # states rather than pinning them (test-dsem-retrospective.R); jitter()
  # perturbs the DSEM starts while leaving the fixed covariate columns alone
  # (test-dsem-jitter.R); self_test() refits the DSEM on simulated observations
  # (test-dsem-self-test.R).
  # remove_F(), reweight_comps() and osa_residuals() are no longer here either.
  # All three refit or read through paths that carry the DSEM in `inits`, which
  # only became true once fit_mod() stopped dropping the dsem_* blocks out of a
  # warm start -- before that they silently rebuilt the model with the
  # recruitment SD at its start value.
  #
  # profile() is no longer a blanket refusal either: a DSEM blocks only the
  # slots it maps out (R_log_sd, rec_dev). See test-dsem-profile.R.
  #
  # process_residuals() still refuses, and for a statistical reason rather than
  # a plumbing one: it standardizes a posterior draw by a per-year normal prior,
  # which is the wrong prior for a GMRF unless the sem is IID.
  testthat::expect_error(Rceattle::process_residuals(fake),
                         "does not support a DSEM")
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
