# Convergence diagnostics (R/0-convergence.R). Exercised on minimal synthetic
# fit-like objects so the suite stays fast (no TMB fits). The real ATF/NORK
# (FAIL) vs pollock (OK) fixtures are the integration-level counterpart.

# Minimal object carrying just what the checks read.
make_fake_fit <- function(R_sd, estimated, max_gradient = 1e-5,
                          worst = "log_F", pdHess = TRUE) {
  structure(list(
    data_list = list(nspp = 1L, spnames = "spp1",
                     random_rec = estimated),
    quantities = list(R_sd = R_sd),
    sdrep = NULL,
    # rec_sd_idx + beta_z map => .rceattle_rec_sd() reads `Estimated` from the
    # map (NA = fixed/off, non-NA = estimated).
    dsem = list(tmb_inputs = list(data = list(rec_sd_idx = 1L))),
    map  = list(mapList = list(beta_z = if (estimated) 1 else NA)),
    .conv_hindcast = list(
      max_gradient = max_gradient,
      worst = list(param = worst, gradient = max_gradient),
      pdHess = pdHess
    )
  ), class = "Rceattle")
}

test_that("estimated recruitment-SD collapse is a FAIL", {
  fit <- make_fake_fit(R_sd = 5e-14, estimated = TRUE,
                       max_gradient = 4e12, worst = "beta_z", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_s3_class(cv, "Rceattle_convergence")
  expect_equal(cv$status, "FAIL")
  expect_true("rec_sd_collapse" %in% names(cv$checks))
  expect_equal(cv$checks$rec_sd_collapse$severity, "FAIL")
  # cross-reference: beta_z worst gradient names the 1/R_sd^2 mechanism
  expect_match(cv$checks$rec_sd_collapse$message, "beta_z")
})

test_that("healthy estimated recruitment SD is OK", {
  fit <- make_fake_fit(R_sd = 0.67, estimated = TRUE,
                       max_gradient = 1e-4, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$status, "OK")
  expect_false("rec_sd_collapse" %in% names(cv$checks))
})

test_that("a fixed (non-estimated) small SD is NOT flagged", {
  # random_rec = FALSE pins sigmaR at the prior; small but legitimate.
  fit <- make_fake_fit(R_sd = 5e-14, estimated = FALSE,
                       max_gradient = 1e-4, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_false("rec_sd_collapse" %in% names(cv$checks))
  expect_equal(cv$status, "OK")
})

test_that("high gradient and non-PD Hessian are flagged", {
  fit <- make_fake_fit(R_sd = 0.5, estimated = TRUE,
                       max_gradient = 4e12, worst = "beta_z", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$status, "FAIL")
  expect_equal(cv$checks$max_gradient$severity, "FAIL")
  expect_equal(cv$checks$pdHess$severity, "FAIL")
})

test_that("Hessian eigen check flags an ill-conditioned covariance", {
  # cov.fixed with one near-flat direction: param p2 ~ 1e8x more uncertain.
  cov <- diag(c(1, 1e8))
  cov[1, 2] <- cov[2, 1] <- 1e3        # correlate so the flat dir mixes p1,p2
  dimnames(cov) <- list(c("beta_z", "beta_z"), c("beta_z", "beta_z"))
  fit <- make_fake_fit(R_sd = 0.5, estimated = TRUE)
  fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_true("hessian_conditioning" %in% names(cv$checks))# fit-like objects so the suite stays fast (no TMB fits).
  expect_true(cv$checks$hessian_conditioning$severity %in% c("WARN", "FAIL"))
  expect_gt(cv$checks$hessian_conditioning$data$condition_number, 1e6)
})


make_fake_fit <- function(max_gradient = 1e-5, worst = "log_F", pdHess = TRUE) {
  structure(list(
    sdrep = NULL,
    .conv_hindcast = list(
      max_gradient = max_gradient,
      worst = list(param = worst, gradient = max_gradient),
      pdHess = pdHess)
  ), class = "Rceattle")
}

test_that("high gradient and non-PD Hessian are flagged", {
  fit <- make_fake_fit(max_gradient = 4e12, worst = "beta_z", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_s3_class(cv, "Rceattle_convergence")
  expect_equal(cv$status, "FAIL")
  expect_equal(cv$checks$max_gradient$severity, "FAIL")
  expect_equal(cv$checks$pdHess$severity, "FAIL")
  expect_match(cv$checks$max_gradient$message, "beta_z")
})

test_that("a converged fit is OK", {
  fit <- make_fake_fit(max_gradient = 1e-5, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$status, "OK")
})

test_that("Hessian eigen check flags an ill-conditioned covariance", {
  cov <- diag(c(1, 1e8))
  cov[1, 2] <- cov[2, 1] <- 1e3
  dimnames(cov) <- list(c("a", "b"), c("a", "b"))
  fit <- make_fake_fit()
  fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_true(cv$checks$hessian_conditioning$severity %in% c("WARN", "FAIL"))
  expect_gt(cv$checks$hessian_conditioning$data$condition_number, 1e6)
})

test_that("Hessian eigen check is OK on a well-conditioned covariance", {
  cov <- diag(c(1, 2))
  dimnames(cov) <- list(c("a", "b"), c("a", "b"))
  fit <- make_fake_fit(R_sd = 0.5, estimated = TRUE)
  fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  # record is present (reports the condition number) but severity is OK
  expect_equal(cv$checks$hessian_conditioning$severity, "OK")
  expect_equal(cv$status, "OK")
})


test_that("pre-fit screen flags sparse + collinear + mis-scaled predictors", {
  set.seed(1)
  yrs <- 1980:2020
  base <- rnorm(length(yrs))
  ed <- data.frame(
    Year = yrs,
    # sst1/sst2 highly collinear; observed only in the last 10 years (sparse)
    sst1 = NA_real_, sst2 = NA_real_,
    # upwelling: fully observed but ~30x larger SD (scale heterogeneity)
    up   = base * 30 + 5
  )
  tail_idx <- which(yrs >= 2011)
  ed$sst1[tail_idx] <- base[tail_idx] + rnorm(length(tail_idx), 0, 0.05)
  ed$sst2[tail_idx] <- ed$sst1[tail_idx] + rnorm(length(tail_idx), 0, 0.05)

  sem_full <- data.frame(
    first  = c("sst1", "sst2", "up"),
    second = c("recdevs1", "recdevs1", "recdevs1"),
    lag    = c(1, 1, 1),
    direction = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  dsem <- list(sem_full = sem_full)
  dl <- list(env_data = ed, styr = 1980, endyr = 2020)

  cv <- check_dsem_spec(dl, dsem)
  expect_s3_class(cv, "Rceattle_convergence")
  expect_true("rec_predictor_observability" %in% names(cv$checks))
  expect_true("covariate_scale" %in% names(cv$checks))
  # sst1/sst2 collinear among the years both are observed
  expect_true("rec_design_conditioning" %in% names(cv$checks))
})

test_that("pre-fit screen is OK for a well-behaved spec", {
  set.seed(2)
  yrs <- 1980:2020
  ed <- data.frame(Year = yrs,
                   a = rnorm(length(yrs)),
                   b = rnorm(length(yrs)))   # independent, fully observed, ~unit SD
  sem_full <- data.frame(first = c("a", "b"),
                         second = c("recdevs1", "recdevs1"),
                         lag = c(1, 1), direction = c(1, 1),
                         stringsAsFactors = FALSE)
  cv <- check_dsem_spec(list(env_data = ed, styr = 1980, endyr = 2020),
                        list(sem_full = sem_full))
  expect_equal(cv$status, "OK")
})

test_that("sdreport failure is a FAIL", {
  fit <- structure(list(.conv_hindcast = list(
    sd_requested = TRUE, sd_present = FALSE)), class = "Rceattle")
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$sdreport_failed$severity, "FAIL")
  # not flagged when sdreport was not requested
  fit2 <- structure(list(.conv_hindcast = list(
    sd_requested = FALSE, sd_present = FALSE)), class = "Rceattle")
  expect_false("sdreport_failed" %in% names(convergence_diagnostics(fit2)$checks))
})

test_that("a parameter at a configured bound is flagged (WARN)", {
  fit <- structure(list(.conv_hindcast = list(
    par   = c(a = 0.5, b = 2.0),     # b sits at its upper bound
    lower = c(0, -1), upper = c(1, 2))), class = "Rceattle")
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$parameters_on_bounds$severity, "WARN")
  expect_match(cv$checks$parameters_on_bounds$message, "b")
  fit2 <- structure(list(.conv_hindcast = list(
    par = c(log_F = -999), lower = c(-999), upper = c(10))), class = "Rceattle")
  expect_false("parameters_on_bounds" %in%
                 names(convergence_diagnostics(fit2)$checks))
})

test_that(".capture_opt_convergence aligns each MLE with its own bounds", {
  # Regression for the bug where the (mle, lower, upper) triple was assembled
  # from two independently-ordered sources, pairing one parameter's MLE with
  # another's bounds: rec_pars (unbounded, MLE 14.9) was reported "at" log_F's
  # [-1000, 10] upper bound. The capture must read the MLE back through the SAME
  # parameter-list shape as the bounds so rows can't drift apart.
  fake_obj <- list(
    env = list(
      last.par.best = c(rec_pars = 14.9, log_F = -999, catch_log_sd = 3),
      random = integer(0),
      # parList(x = par[lfixed()], par = last.par) expands the reduced vector
      # into the full param-list shape that `bounds` and `mapFactor` share.
      parList = function(x = NULL, par = NULL) list(rec_pars = 14.9,
                                                    log_F = -999,
                                                    catch_log_sd = 3)
    ),
    gr = function(p) rep(0, length(p))
  )
  bounds <- list(
    lower = list(rec_pars = -Inf, log_F = -1000, catch_log_sd = -10),
    upper = list(rec_pars =  Inf, log_F =    10, catch_log_sd =   3)
  )
  mapFactor <- list(rec_pars = 1, log_F = 1, catch_log_sd = 1)

  snap <- .capture_opt_convergence(opt = list(), obj = fake_obj,
                                   bounds = bounds, mapFactor = mapFactor,
                                   random_vars = NULL, getsd = FALSE)
  # bounds line up with their own parameter
  expect_equal(unname(snap$par["rec_pars"]), 14.9)
  expect_equal(snap$lower[match("rec_pars", names(snap$par))], -Inf)
  expect_equal(snap$upper[match("catch_log_sd", names(snap$par))], 3)

  fit <- structure(list(.conv_hindcast = snap), class = "Rceattle")
  cv  <- convergence_diagnostics(fit)
  # rec_pars is unbounded -> must NOT be flagged; log_F is a -999 sentinel ->
  # skipped; only catch_log_sd genuinely sits at its upper bound.
  ob <- cv$checks$parameters_on_bounds
  expect_equal(ob$severity, "WARN")
  expect_match(ob$message, "catch_log_sd")
  expect_false(grepl("rec_pars", ob$message))
  expect_false(grepl("log_F", ob$message))
})

test_that("a phase that ends with a high gradient is flagged (WARN)", {
  fit <- structure(list(.conv_phase = list(
    list(phase = 1, max_grad = 1e-4),
    list(phase = 2, max_grad = 5e3))), class = "Rceattle")
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$phasing$severity, "WARN")
  expect_match(cv$checks$phasing$message, "phase 2")
})

test_that("print method runs and is non-erroring", {
  fit <- make_fake_fit(max_gradient = 4e12, worst = "beta_z", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_output(print(cv), "status: FAIL")
  expect_invisible(print(cv))
})
