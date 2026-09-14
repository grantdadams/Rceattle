# A linkage intercept fixed at its init (est_phase = 0) must hold that value when
# fit_mod() is given `inits`. Before 5.33.0 the init was applied only to fresh
# parameters, so a refit at 500 warm-started from a fit at 1170 stayed at 1170.

.ws_fc <- function() fit_control(getsd = FALSE, verbose = 0)

test_that("a fixed alpha linkage holds its init over supplied inits", {
  fixed_at <- function(v) build_srr(srr_fun = "BevertonHolt", linkages = list(
    alpha = linkage_spec(~ 1, est_phase = 0, init = list(`(Intercept)` = v))))
  m1 <- suppressWarnings(fit_mod(Rceattle::BS2017SS, recFun = fixed_at(1170),
                                 estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  m2 <- suppressWarnings(fit_mod(Rceattle::BS2017SS, inits = m1$estimated_params,
                                 recFun = fixed_at(500),
                                 estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  expect_true(all(is.na(m2$map$mapList$rec_pars[, 2])))
  expect_equal(unname(exp(m2$obj$env$parList()$rec_pars[, 2])), rep(500, 3))
})

test_that("an estimable intercept still starts from the supplied inits", {
  est_at <- function(v) build_srr(srr_fun = "BevertonHolt", linkages = list(
    alpha = linkage_spec(~ 1, init = list(`(Intercept)` = v))))
  m1 <- suppressWarnings(fit_mod(Rceattle::BS2017SS, recFun = est_at(1170),
                                 estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  m2 <- suppressWarnings(fit_mod(Rceattle::BS2017SS, inits = m1$estimated_params,
                                 recFun = est_at(500),
                                 estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  expect_equal(unname(exp(m2$obj$env$parList()$rec_pars[, 2])), rep(1170, 3))
})

test_that("a fixed alpha linkage also wins over srr_alpha_init", {
  rf <- build_srr(srr_fun = "BevertonHolt", srr_alpha_init = 50, linkages = list(
    alpha = linkage_spec(~ 1, est_phase = 0, init = list(`(Intercept)` = 1170))))
  m <- suppressWarnings(fit_mod(Rceattle::BS2017SS, recFun = rf,
                                estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  expect_equal(unname(exp(m$obj$env$parList()$rec_pars[, 2])), rep(1170, 3))
})

test_that("a refit keeps the fitted alpha, not the starting-value overrides", {
  rf <- build_srr(srr_fun = "BevertonHolt", srr_alpha_init = 50, linkages = list(
    alpha = linkage_spec(~ 1, est_phase = 0, init = list(`(Intercept)` = 1170))))
  m <- suppressWarnings(fit_mod(Rceattle::BS2017SS, recFun = rf,
                                estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  r <- suppressWarnings(Rceattle:::.refit_like(data_list = m$data_list,
                                               inits = m$estimated_params, estimateMode = 3))
  expect_equal(unname(exp(r$obj$env$parList()$rec_pars[, 2])), rep(1170, 3))
})

test_that("srr_est_mode = 'Fixed' with a fixed alpha linkage is refused", {
  expect_error(suppressWarnings(build_srr(
    srr_fun = "BevertonHolt", srr_est_mode = "Fixed", srr_prior = 20,
    linkages = list(alpha = linkage_spec(~ 1, est_phase = 0, init = list(`(Intercept)` = 30))))),
    "both fix alpha")
})

test_that("profile() of a linkage-fixed M moves species 1's M to each grid value", {
  testthat::skip_on_cran()
  fit <- suppressWarnings(suppressMessages(fit_mod(
    Rceattle::BS2017SS, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    M1Fun = build_M1(M1_model = 1, linkages = list(
      M1 = linkage_spec(~ 1, est_phase = 0, init = list(`(Intercept)` = 0.3)))),
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  pr <- suppressWarnings(suppressMessages(stats::profile(
    fit, param = "log_M1", slots = list(c(1, 1, 1)), values = list(c(0.2, 0.4)))))
  nll <- pr$nll
  expect_length(nll, 2L)
  expect_gt(abs(diff(nll)), 1)
})

test_that("a fixed M linkage holds its init over supplied inits", {
  fixed_M <- function(v) build_M1(M1_model = 1, linkages = list(
    M1 = linkage_spec(~ 1, est_phase = 0, init = list(`(Intercept)` = v))))
  m1 <- suppressWarnings(fit_mod(Rceattle::BS2017SS, M1Fun = fixed_M(0.3),
                                 estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  m2 <- suppressWarnings(fit_mod(Rceattle::BS2017SS, inits = m1$estimated_params,
                                 M1Fun = fixed_M(0.2),
                                 estimateMode = 3, msmMode = 0, fit_control = .ws_fc()))
  expect_equal(unique(round(as.numeric(exp(m2$obj$env$parList()$log_M1[, 1, 1])), 12)), 0.2)
})
