# Ianelli configuration (srr_fun < 2 with srr_pred_fun > 1): the hindcast is
# mean recruitment and the stock-recruit curve enters as a lognormal penalty.
# The curve's alpha is still estimated, so its steepness must be derived and
# checked -- keying either on srr_fun alone leaves the penalty curve unexamined.

test_that("steepness is derived under the Ianelli configuration", {
  m <- fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, msmMode = 0,
    estimateMode = 3,
    recFun = build_srr(srr_fun = 0, srr_pred_fun = 2, srr_est_mode = 1),
    fit_control = fit_control(getsd = FALSE, verbose = 0))

  h <- m$quantities$steepness[, 1]
  # A placeholder (the pre-fix behaviour) sat at 0.99 for every species.
  expect_false(all(abs(h - 0.99) < 1e-8))

  # h = alpha * SPR0 / (4 + alpha * SPR0). alpha is not REPORTed, so read it
  # off the parameter it is stored on: rec_pars[, 2] holds log(alpha).
  alpha <- exp(m$estimated_params$rec_pars[, 2])
  spr0  <- as.numeric(m$quantities$SPR0)
  expect_equal(as.numeric(h), as.numeric(alpha * spr0 / (4 + alpha * spr0)),
               tolerance = 1e-8)
})

test_that("a sub-replacement alpha is flagged under the Ianelli configuration", {
  # srr_prior is an ALPHA when srr_est_mode = 0. 0.8 is below 1/SPR0 for two of
  # the three BS2017SS species, so the curve cannot replace itself.
  m <- suppressWarnings(fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, msmMode = 0,
    estimateMode = 3,
    recFun = build_srr(srr_fun = 0, srr_pred_fun = 2, srr_est_mode = 0,
                       srr_prior = 0.8),
    fit_control = fit_control(getsd = FALSE, verbose = 0)))

  checks <- convergence_diagnostics(m)$checks
  # Gated on srr_fun this check returned an empty list and the fit passed silently.
  expect_true("stock_recruit" %in% names(checks))
  expect_match(checks$stock_recruit$message, "steepness below 0\\.2")
  expect_true(any(m$quantities$steepness[, 1] < 0.2))
})

test_that("a well-posed Ianelli curve reports no stock-recruit problem", {
  m <- fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, msmMode = 0,
    estimateMode = 3,
    recFun = build_srr(srr_fun = 0, srr_pred_fun = 2, srr_est_mode = 1),
    fit_control = fit_control(getsd = FALSE, verbose = 0))

  sr <- convergence_diagnostics(m)$checks$stock_recruit
  expect_true(!is.null(sr))
  expect_match(sr$message, "well posed")
})

test_that("mean recruitment with no penalty curve is still skipped", {
  m <- fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, msmMode = 0,
    estimateMode = 3,
    recFun = build_srr(srr_fun = 0, srr_pred_fun = 0),
    fit_control = fit_control(getsd = FALSE, verbose = 0))

  expect_false("stock_recruit" %in% names(convergence_diagnostics(m)$checks))
})

test_that("build_srr warns when a steepness is passed as srr_est_mode 0's alpha", {
  expect_warning(
    build_srr(srr_fun = 0, srr_pred_fun = 2, srr_est_mode = 0, srr_prior = 0.8),
    "is an alpha here, not a steepness")

  # A plausible alpha must not warn, and modes 2/3 genuinely do take a steepness.
  expect_no_warning(
    build_srr(srr_fun = 0, srr_pred_fun = 2, srr_est_mode = 0, srr_prior = 1170))
  expect_no_warning(
    build_srr(srr_fun = 0, srr_pred_fun = 2, srr_est_mode = 2,
              srr_prior = 0.8, srr_prior_sd = 0.2))
  # Ricker's srr_prior is an alpha in every mode, so there is nothing to confuse.
  expect_no_warning(
    build_srr(srr_fun = 0, srr_pred_fun = 4, srr_est_mode = 0, srr_prior = 0.8))
})
