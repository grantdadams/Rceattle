# Convergence diagnostics (R/0-convergence.R). Exercised on minimal synthetic
# fit-like objects so the suite stays fast (no TMB fits).

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
  fit <- make_fake_fit(max_gradient = 4e12, worst = "sel_inf", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_s3_class(cv, "Rceattle_convergence")
  expect_equal(cv$status, "FAIL")
  expect_equal(cv$checks$max_gradient$severity, "FAIL")
  expect_equal(cv$checks$pdHess$severity, "FAIL")
  expect_match(cv$checks$max_gradient$message, "sel_inf")
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

test_that("Hessian eigen check names a diffusely-loaded flat direction", {
  # Flat direction spread evenly over 40 rec_dev coefficients (each loading
  # ~0.16, below any single-coefficient threshold) plus a little ln_srv_sel.
  # The pre-fix check printed "loads on: ." with nothing after it.
  nm  <- c(rep("rec_dev", 40), rep("ln_srv_sel", 4))
  p   <- length(nm)
  v1  <- c(rep(1, 40), rep(1.5, 4)); v1 <- v1 / sqrt(sum(v1^2))
  Q   <- qr.Q(qr(cbind(v1, diag(p)[, -1])))
  cov <- Q %*% diag(c(1e8, rep(1, p - 1))) %*% t(Q)
  dimnames(cov) <- list(nm, nm)
  fit <- make_fake_fit(); fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  hc  <- convergence_diagnostics(fit)$checks$hessian_conditioning
  expect_true(hc$severity %in% c("WARN", "FAIL"))
  expect_match(hc$message, "loads on: [A-Za-z]")   # never blank after "loads on: "
  expect_match(hc$message, "rec_dev")              # the dominant block is named
  expect_true("rec_dev" %in% hc$data$loadings$param)
})

test_that("Hessian eigen check falls back to par.fixed names without dimnames", {
  nm  <- c("a", "b", "c")
  cov <- diag(c(1e8, 1, 1))                        # flat direction is coeff 'a'
  fit <- make_fake_fit()
  fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE,
                    par.fixed = stats::setNames(c(0, 0, 0), nm))
  hc  <- convergence_diagnostics(fit)$checks$hessian_conditioning
  expect_match(hc$message, "loads on: a")          # named from par.fixed, not "p1"
})

test_that("Hessian eigen check is OK on a well-conditioned covariance", {
  cov <- diag(c(1, 2)); dimnames(cov) <- list(c("a", "b"), c("a", "b"))
  fit <- make_fake_fit(); fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$hessian_conditioning$severity, "OK")
  expect_equal(cv$status, "OK")
})

test_that("sdreport failure is a FAIL", {
  fit <- structure(list(.conv_hindcast = list(
    sd_requested = TRUE, sd_present = FALSE)), class = "Rceattle")
  expect_equal(convergence_diagnostics(fit)$checks$sdreport_failed$severity,
               "FAIL")
  fit2 <- structure(list(.conv_hindcast = list(
    sd_requested = FALSE, sd_present = FALSE)), class = "Rceattle")
  expect_false("sdreport_failed" %in%
                 names(convergence_diagnostics(fit2)$checks))
})

test_that("a parameter at a configured bound is flagged (WARN)", {
  fit <- structure(list(.conv_hindcast = list(
    par   = c(a = 0.5, b = 2.0),
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
  fit <- make_fake_fit(max_gradient = 4e12, worst = "sel_inf", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_output(print(cv), "status: FAIL")
  expect_invisible(print(cv))
})
