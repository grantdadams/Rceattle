# IID random-effect linkage density: `~ (1 | group)` expands to one deviation
# per level, damped by a N(0, sigma) density in jnll_comp row 20 (0-based; the
# "Linkage random effects" row). These cover the registry/encoding wiring
# (fast) and the end-to-end fit + density (skip_on_cran).

# ---- registry + encoding (fast, no fit) -------------------------------------

testthat::test_that("an IID linkage builds a bijective RE registry", {
  env <- data.frame(Year = 2000:2009)
  spec <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet,
                                 fleet = 1L, link = "log")
  pool <- Rceattle:::pool_linkages(list(q = list(q = spec)), env_data = env,
                                   strata = list(fleet = 1L))
  tbl <- pool$table
  re  <- tbl[!is.na(tbl$re_index), ]

  testthat::expect_equal(nrow(re), 10L)                 # one slot per year
  testthat::expect_setequal(re$re_index, 0:9)           # bijection
  testthat::expect_true(all(re$sigma_index == 0L))      # single group
  testthat::expect_false(anyNA(re$re_time))             # numeric time captured
  # fixed rows carry no registry entry
  testthat::expect_true(all(is.na(tbl$re_index[is.na(tbl$re_struct)])))
})

testthat::test_that("distinct fleets get distinct sigma groups", {
  env <- data.frame(Year = 2000:2004)
  s1 <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet, fleet = 1L)
  s2 <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet, fleet = 2L)
  pool <- Rceattle:::pool_linkages(list(q = list(q = list(s1, s2))),
                                   env_data = env, strata = list(fleet = 1:2))
  enc <- Rceattle:::encode_linkage_for_tmb(pool$table, pool$X)
  re  <- pool$table[!is.na(pool$table$re_index), ]

  testthat::expect_equal(max(re$sigma_index) + 1L, 2L)  # two groups
  # per-slot sigma group matches the owning row
  expect_sigma <- integer(nrow(re))
  expect_sigma[re$re_index + 1L] <- re$sigma_index
  testthat::expect_identical(enc$linkage_re_sigma, as.integer(expect_sigma))
})

testthat::test_that("a fixed-effect prior on an RE column is rejected", {
  env <- data.frame(Year = 2000:2004)
  spec <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet, fleet = 1L)
  pool <- Rceattle:::pool_linkages(list(q = list(q = spec)), env_data = env,
                                   strata = list(fleet = 1L))
  bad <- pool$table
  i <- which(!is.na(bad$re_struct))[1]
  bad$prior_family[i] <- "normal"; bad$prior_p1[i] <- 0; bad$prior_p2[i] <- 1
  testthat::expect_error(
    Rceattle:::encode_linkage_for_tmb(bad, pool$X), "fixed-effect prior")
})

# ---- end-to-end density (skip_on_cran) --------------------------------------

.re_row <- function(fit) {
  which(rownames(fit$quantities$jnll_comp) == "Linkage random effects")
}

testthat::test_that("the RE density equals -sum(dnorm(dev, 0, sigma)) exactly", {
  # The strongest check of the density math, and it needs no optimisation:
  # build the object (estimateMode = 3), evaluate it once at a KNOWN deviate
  # pattern + SD, and confirm the reported density row is the mathematical
  # definition to machine precision. (The legacy M1_re = 2 density is the same
  # SCALE(AR1(rho = 0), sigma) form, i.e. the identical dnorm sum.)
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7L)))
  obj <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
    qFun = qfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))$obj

  p <- obj$env$last.par
  n_re <- sum(names(p) == "beta_linkage_re")
  set.seed(1)
  dev  <- stats::rnorm(n_re, 0, 0.4)
  lsig <- log(0.3)
  p[names(p) == "beta_linkage_re"]   <- dev
  p[names(p) == "log_sigma_linkage"] <- lsig
  rep <- obj$report(p)
  # raw REPORT has no rownames; the density is the last jnll_comp row
  # (cpp row 20, 0-based).
  reported <- sum(rep$jnll_comp[nrow(rep$jnll_comp), ])
  testthat::expect_equal(reported,
    -sum(stats::dnorm(dev, 0, exp(lsig), log = TRUE)), tolerance = 1e-8)
})

testthat::test_that("a model without a random linkage keeps density row 20 at 0", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  testthat::expect_equal(sum(fit$quantities$jnll_comp[.re_row(fit), ]), 0)
})

testthat::test_that("an IID (1|Year) q linkage builds finite and estimates a density", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS                      # env_data already has Year
  survey_flt <- 7L  # EIT_Pollock, the only Estimated-q survey (Fixed q is rejected)
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = survey_flt)))
  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)

  # build-only: obj$fn / obj$gr finite (no NaN in the Laplace tape)
  b3 <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                          qFun = qfun, fit_control = ctl)
  testthat::expect_true(is.finite(b3$obj$fn()))
  testthat::expect_true(all(is.finite(b3$obj$gr())))

  # fit: one deviation per hindcast year, an active finite density, finite sigma
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                           qFun = qfun, fit_control = ctl)
  n_hind <- d$endyr - d$styr + 1L
  testthat::expect_equal(length(fit$estimated_params$beta_linkage_re), n_hind)
  testthat::expect_true(sum(fit$quantities$jnll_comp[.re_row(fit), ]) != 0)
  sigma <- exp(fit$estimated_params$log_sigma_linkage)
  testthat::expect_true(is.finite(sigma) && sigma > 0)
})

testthat::test_that("rw(1|Year) fits, pins the first deviate, and uses a difference density", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 7L)))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  dev   <- fit$estimated_params$beta_linkage_re
  sigma <- exp(fit$estimated_params$log_sigma_linkage)
  # first deviate pinned at 0 (the walk's level is carried by the base q)
  testthat::expect_equal(dev[1], 0, tolerance = 1e-10)
  # the density is N(0, sigma) on successive first differences (deviates are in
  # elapsed-time order)
  testthat::expect_equal(sum(fit$quantities$jnll_comp[.re_row(fit), ]),
    -sum(stats::dnorm(diff(dev), 0, sigma, log = TRUE)), tolerance = 1e-7)
})

testthat::test_that("ar1(1|Year) fits with a stationary AR1 density and estimable rho", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 7L)))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  re    <- fit$estimated_params$beta_linkage_re
  sigma <- exp(fit$estimated_params$log_sigma_linkage)
  rho   <- tanh(fit$estimated_params$trans_rho_linkage)   # rho_trans(x) == tanh(x)
  testthat::expect_true(is.finite(rho) && abs(rho) < 1)
  # exact stationary AR1 negative log density, MARGINAL SD = sigma:
  n <- length(re)
  nll <- -stats::dnorm(re[1], 0, sigma, log = TRUE) -
          sum(stats::dnorm(re[2:n], rho * re[1:(n - 1)], sigma * sqrt(1 - rho^2), log = TRUE))
  testthat::expect_equal(sum(fit$quantities$jnll_comp[.re_row(fit), ]), nll, tolerance = 1e-6)
})

testthat::test_that("init = list(rho = v) fixes the ar1 correlation at v", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 7L,
                               init = list(rho = 0.5))))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  testthat::expect_equal(tanh(fit$estimated_params$trans_rho_linkage), 0.5, tolerance = 1e-6)
})

testthat::test_that("rw()/ar1() require a numeric grouping variable", {
  env <- data.frame(Year = letters[1:5])
  testthat::expect_error(
    Rceattle:::pool_linkages(list(q = list(q = Rceattle::linkage_spec(
      ~ rw(1 | Year), param = "q", by = ~ fleet, fleet = 1L))),
      env_data = env, strata = list(fleet = 1L)),
    "numeric grouping")
})

testthat::test_that("init = list(sigma = v) fixes the deviation SD at v", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  survey_flt <- 7L  # EIT_Pollock, the only Estimated-q survey (Fixed q is rejected)
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = survey_flt,
                               init = list(sigma = 0.1))))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  # the SD is held exactly at the input value, and the density is still active
  testthat::expect_equal(exp(fit$estimated_params$log_sigma_linkage), 0.1, tolerance = 1e-9)
  testthat::expect_true(sum(fit$quantities$jnll_comp[.re_row(fit), ]) != 0)
})

testthat::test_that("priors = list(sigma = ...) puts a prior on the deviation SD", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  survey_flt <- 7L  # EIT_Pollock, the only Estimated-q survey (Fixed q is rejected)
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = survey_flt,
                               priors = list(sigma = lognormal(log(0.05), 0.3)))))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  pr_row <- which(rownames(fit$quantities$jnll_comp) == "Linkage-table priors")
  testthat::expect_true(sum(fit$quantities$jnll_comp[pr_row, ]) != 0)   # prior contributes
  testthat::expect_true(is.finite(exp(fit$estimated_params$log_sigma_linkage)))
})
