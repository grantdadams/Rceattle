# Estimate_catch_sd = "Analytical" (2) is the Ludwig and Walters (1994)
# concentrated estimator: sigma_hat = sqrt( sum(log(obs) - log(pred))^2 / n ).
#
# Until 5.12.0 the option was documented in the schema, accepted by
# validate_switches(), and treated as a real mode by build_map() (which maps
# catch_log_sd out) -- but the template's dispatch had only case 0 and case 1,
# so the fit died with "Invalid 'Estimate_sigma_catch'" after passing every R
# check. est_sigma_index had implemented case 2 all along; this is its mirror.

test_that("the analytical catch sd equals the closed form", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  d <- BS2017SS
  # NOTE: Fleet_type is integer-coded on the raw bundled data (1 = Fishery).
  # Subsetting on the string "Fishery" here matches nothing and silently
  # produces a mode-0 model -- which is how a first draft of this test passed
  # while exercising nothing.
  d$fleet_control$Estimate_catch_sd[d$fleet_control$Fleet_type == 1] <- 2

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))))

  # The switch really reached the template.
  testthat::expect_true(all(fit$obj$env$data$est_sigma_fsh[1:3] == 2))

  dl <- fit$obj$env$data
  cc <- dl$catch_ctl; co <- dl$catch_obs
  hat <- fit$quantities$catch_hat
  got <- fit$quantities$log_catch_analytical_sd

  checked <- 0L
  for (flt in sort(unique(cc[, 1]))) {
    if (dl$flt_type[flt] != 1) next
    rows <- which(cc[, 1] == flt & cc[, 3] > 0 & cc[, 3] <= dl$endyr & co[, 1] > 0)
    if (!length(rows)) next
    checked <- checked + 1L
    want <- sqrt(sum((log(co[rows, 1]) - log(hat[rows]))^2) / length(rows))
    testthat::expect_equal(as.numeric(got[flt]), want, tolerance = 1e-10,
                           info = paste("fleet", flt))
  }
  testthat::expect_gt(checked, 0L)
})

test_that("the analytical sd is the one the likelihood actually applies", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  d <- BS2017SS
  d$fleet_control$Estimate_catch_sd[d$fleet_control$Fleet_type == 1] <- 2

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                                          verbose = 0))))
  # Computing it and then not using it would look identical from the outside.
  applied <- fit$quantities$catch_sd
  applied <- sort(unique(round(applied[applied > 0], 8)))
  analytical <- sort(round(fit$quantities$log_catch_analytical_sd[1:3], 8))
  testthat::expect_equal(applied, as.numeric(analytical))

  testthat::expect_true(is.finite(fit$opt$objective))
})

test_that("a model that does not ask for it is bit-identical", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The whole point: a new dispatch arm must not perturb the modes that were
  # already there. BS2017SS ships at Estimate_catch_sd = 0.
  data("BS2017SS")
  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 1, msmMode = 0,
                      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                                          verbose = 0))))
  # The pinned value for this configuration; see test-golden-regression.R for
  # the phased/polished references.
  testthat::expect_equal(fit$opt$objective, 10241.03042750, tolerance = 1e-6)
})
