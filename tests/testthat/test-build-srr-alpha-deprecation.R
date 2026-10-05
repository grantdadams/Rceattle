# Alpha priors, fixed values and starting values outside a linkage are deprecated
# (5.33.0); srr_prior stays for the Beverton-Holt steepness prior.

test_that("each deprecated alpha input warns with its replacement", {
  expect_warning(build_srr(srr_fun = "BevertonHolt", srr_est_mode = "Fixed", srr_prior = 20),
                 "est_phase = 0")
  expect_warning(build_srr(srr_fun = "Ricker", srr_est_mode = "LognormalPrior",
                           srr_prior = 5, srr_prior_sd = 0.5),
                 "prior_lognormal")
  expect_warning(build_srr(srr_fun = "Ricker", srr_prior = 5), "srr_alpha_init")
  expect_warning(build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt", srr_prior = 5),
                 "srr_alpha_init")
})

test_that("the steepness prior, the defaults and the replacements do not warn", {
  expect_no_warning(build_srr())
  expect_no_warning(build_srr(srr_fun = "Ricker"))
  expect_no_warning(build_srr(srr_fun = "BevertonHolt", srr_est_mode = "LognormalPrior",
                              srr_prior = 0.8, srr_prior_sd = 0.2))
  expect_no_warning(build_srr(srr_fun = "BevertonHolt", srr_est_mode = "BetaPrior",
                              srr_prior = 0.8, srr_prior_sd = 0.2))
  expect_no_warning(build_srr(srr_fun = "Ricker", srr_alpha_init = 5))
  expect_no_warning(build_srr(srr_fun = "Ricker", linkages = list(alpha = linkage_spec(
    ~ 1, priors = list(`(Intercept)` = prior_lognormal(log(5), 0.5))))))
})

test_that("a deprecated alpha prior with a linkage prior is still refused", {
  expect_error(suppressWarnings(build_srr(
    srr_fun = "Ricker", srr_est_mode = "LognormalPrior", srr_prior = 5, srr_prior_sd = 0.5,
    linkages = list(alpha = linkage_spec(
      ~ 1, priors = list(`(Intercept)` = prior_lognormal(log(5), 0.5)))))),
    "use one, not both")
})
