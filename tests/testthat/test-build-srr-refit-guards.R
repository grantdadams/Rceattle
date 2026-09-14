# build_srr() / refit guards. A Ricker LognormalPrior and an alpha intercept prior
# are the same density (counted twice, 0.7711 each). The map warning fired on every
# refit's own fitted map; srr_pred_fun 1 under srr_fun 0 refit without its penalty.

.ricker_link_prior <- function(key = "(Intercept)") {
  pr <- list(prior_lognormal(log(5), 0.5)); names(pr) <- key
  list(alpha = linkage_spec(~ 1, priors = pr))
}

test_that("a Ricker LognormalPrior plus an alpha intercept prior is refused", {
  for (key in c("(Intercept)", "intercept")) {
    expect_error(suppressWarnings(build_srr(srr_fun = "Ricker", srr_est_mode = "LognormalPrior",
                                            srr_prior = 5, srr_prior_sd = 0.5,
                                            linkages = .ricker_link_prior(key)),
                                  classes = "rceattle_deprecated"),
                 "use one, not both", info = key)
  }
  # Either one alone is fine, and so is a Beverton-Holt steepness prior with it.
  expect_no_error(suppressWarnings(build_srr(srr_fun = "Ricker", linkages = .ricker_link_prior())))
  expect_no_error(suppressWarnings(build_srr(
    srr_fun = "BevertonHolt", srr_est_mode = "LognormalPrior", srr_prior = 0.8,
    srr_prior_sd = 0.2, linkages = .ricker_link_prior())))
})

test_that("a retired srr_pred_fun = 1 under the penalty form says its penalty is dropped", {
  expect_warning(m <- Rceattle:::.srr_fun_structural(1L, penalty = TRUE),
                 "without that penalty")
  expect_identical(m, 0L)
  expect_warning(Rceattle:::.srr_fun_structural(1L), "no effect since 4.4.0")
})

.map_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  grep("fixes both stock-recruit parameters", w, value = TRUE)
}

test_that("the map warning flags a stuck curve, and only a stuck curve", {
  set.seed(123)
  d <- make_msm_test_data()$data_list
  build <- function(map, recFun) suppressMessages(fit_mod(
    data_list = d, inits = NULL, map = map, estimateMode = 3, msmMode = 1,
    suitMode = 0, initMode = "NonEquilibrium", random_rec = FALSE, recFun = recFun,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))
  pen <- function(alpha = NULL) build_srr(
    srr_fun = "mean", srr_pred_fun = "BevertonHolt",
    linkages = if (!is.null(alpha)) list(alpha = alpha))
  # The model's own map with alpha and beta off. (A mean-recruitment map cannot be
  # carried onto a linkage model at all: fit_mod() stops on the beta_linkage size.)
  stuck <- function(alpha) {
    m <- suppressWarnings(build(NULL, pen(alpha)))$map
    m$mapList$rec_pars[, 2:3] <- NA
    m$mapFactor$rec_pars <- factor(m$mapList$rec_pars)
    m
  }

  # An intercept-only alpha linkage keeps its level on the mapped rec_pars: still stuck.
  a1 <- linkage_spec(~ 1)
  expect_length(.map_warnings(build(stuck(a1), pen(a1))), 1L)

  # An estimated slope moves the curve, for every species it covers.
  a2 <- linkage_spec(~ 0 + EnvData)
  expect_length(.map_warnings(build(stuck(a2), pen(a2))), 0L)
  # Covering species 1 only leaves species 2 stuck.
  a3 <- linkage_spec(~ 0 + EnvData, species = 1)
  w3 <- .map_warnings(build(stuck(a3), pen(a3)))
  expect_match(w3, d$spnames[2], fixed = TRUE)
  expect_no_match(w3, d$spnames[1], fixed = TRUE)

  # A debug or projection map has every rec_pars entry off.
  off <- suppressWarnings(build(NULL, build_srr()))$map
  off$mapList$rec_pars[] <- NA
  off$mapFactor$rec_pars <- factor(off$mapList$rec_pars)
  expect_length(.map_warnings(build(off, pen())), 0L)
})
