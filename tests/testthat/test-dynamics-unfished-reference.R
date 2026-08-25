# Multispecies unfished SSB (`MSSB0`) has to survive on the fitted object.
#
# Under `msmMode > 0` the template discards its own equilibrium SB0 and reads
# the `MSSB0` DATA_VECTOR instead (src/TMB/ceattle.cpp:2105-2109), so
# `ssb_depletion` is `ssb / MSSB0`. `clean_data()` seeds `MSSB0` at a 999 mt
# placeholder because no workbook can supply it, and fit_mod() section 10.2
# derives the real value by projecting under no fishing.
#
# Until 5.15.0 that derived value was written only into the reorganized copy
# section 10.2 refits from, so the *returned* `data_list` kept the 999. Every
# refit off a fitted object -- `.refit_like()`, `remove_F()`, and each
# `run_mse()` projection -- re-entered the template with the placeholder and
# read SSB/999 as depletion. Measured on the Pacific hake three-species model,
# terminal depletion came back as 2.68e3 where the dynamic reference gave 0.96.

test_that("DynamicHCR gives a multispecies model a dynamic-B0 depletion", {
  testthat::skip_on_cran()

  # With HCR = 0 the projection is unfished, so terminal-year SSB is
  # multispecies SB0 and ssb_depletion is measured against it. That arm used to
  # run unconditionally and overwrote the DynamicHCR result, pinning terminal
  # depletion at exactly 1 and leaving no way to get a dynamic B0 out of a
  # multispecies model.
  set.seed(123)
  simData <- make_msm_test_data()$data_list

  fit <- suppressMessages(fit_mod(
    data_list    = simData,
    inits        = NULL,
    estimateMode = 4,
    msmMode      = 1,
    suitMode     = 0,
    initMode     = "NonEquilibrium",
    random_rec   = FALSE,
    HCR          = build_hcr(HCR = 0, DynamicHCR = TRUE),
    fit_control  = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)
  ))

  dep <- fit$quantities$ssb_depletion
  expect_equal(dep, fit$quantities$ssb / fit$quantities$DynamicSB0,
               tolerance = 1e-12)

  # And specifically not the terminal-year index, whose last column is all 1.
  expect_false(isTRUE(all.equal(unname(dep[, ncol(dep)]),
                                rep(1, nrow(dep)), tolerance = 1e-8)))
})


test_that("clean_data() seeds MSSB0 with the placeholder mse_summary() tests against", {
  # Two files have to agree on this value or the summary silently reports
  # SSB/999 as a depletion again: R/0-clean_data.R seeds it, R/10-mse_summary.R
  # reads it back to decide the equilibrium reference was never derived.
  set.seed(123)
  simData <- make_msm_test_data()$data_list
  simData$MSSB0 <- NULL
  simData$MSB0 <- NULL

  cleaned <- suppressMessages(Rceattle:::clean_data(simData))

  expect_equal(cleaned$MSSB0, rep(Rceattle:::.RCE_MSSB0_PLACEHOLDER, cleaned$nspp))
  expect_equal(cleaned$MSB0,  rep(Rceattle:::.RCE_MSSB0_PLACEHOLDER, cleaned$nspp))
})


test_that("a multispecies HCR fit carries its derived MSSB0 onto the returned data_list", {
  testthat::skip_on_cran()

  set.seed(123)
  simData <- make_msm_test_data()$data_list

  fit <- suppressMessages(fit_mod(
    data_list    = simData,
    inits        = NULL,
    estimateMode = 4,             # projection only; section 10.2 still derives MSSB0
    msmMode      = 1,
    suitMode     = 0,
    initMode     = "NonEquilibrium",
    random_rec   = FALSE,
    HCR          = build_hcr(HCR = 6, Ptarget = 0.4, Plimit = 0.2),
    fit_control  = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)
  ))

  # The placeholder must not survive the fit.
  expect_false(any(fit$data_list$MSSB0 == 999),
               info = "MSSB0 is still the 999 mt placeholder on the returned data_list")
  expect_false(any(fit$data_list$MSB0 == 999),
               info = "MSB0 is still the 999 mt placeholder on the returned data_list")

  # And it must be the value the template actually used: under msmMode > 0 the
  # template sets SB0(sp, yr) = MSSB0(sp) for every year, so the terminal
  # column of the reported SB0 is MSSB0 echoed back.
  sb0 <- fit$quantities$SB0
  expect_equal(as.numeric(fit$data_list$MSSB0),
               as.numeric(sb0[, ncol(sb0)]),
               tolerance = 1e-10)
})


test_that("a refit from a fitted multispecies object reproduces its SB0", {
  testthat::skip_on_cran()

  set.seed(123)
  simData <- make_msm_test_data()$data_list

  fit <- suppressMessages(fit_mod(
    data_list    = simData,
    inits        = NULL,
    estimateMode = 4,
    msmMode      = 1,
    suitMode     = 0,
    initMode     = "NonEquilibrium",
    random_rec   = FALSE,
    HCR          = build_hcr(HCR = 6, Ptarget = 0.4, Plimit = 0.2),
    fit_control  = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)
  ))

  # Refitting from the object's own data list is the contract `.refit_like()`,
  # `remove_F()` and run_mse()'s projections all rely on: the same model in,
  # the same unfished reference out. With the placeholder left in place this
  # returned 999 regardless of what the first fit derived.
  refit <- suppressMessages(fit_mod(
    data_list    = fit$data_list,
    inits        = fit$estimated_params,
    estimateMode = 4,
    msmMode      = 1,
    suitMode     = 0,
    initMode     = "NonEquilibrium",
    random_rec   = FALSE,
    HCR          = build_hcr(HCR = 6, Ptarget = 0.4, Plimit = 0.2),
    fit_control  = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)
  ))

  # Assert the shared value is a real one first: two fits that both lost it
  # agree on 999 trivially, so equality alone would pass with the bug present.
  expect_false(any(refit$data_list$MSSB0 == Rceattle:::.RCE_MSSB0_PLACEHOLDER),
               info = "refit fell back to the 999 mt placeholder")
  expect_equal(as.numeric(refit$data_list$MSSB0),
               as.numeric(fit$data_list$MSSB0),
               tolerance = 1e-10)
})
