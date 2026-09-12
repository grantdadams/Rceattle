# data_check() stock-recruit guards. Penalty years past endyr scored projected
# recruitment (BS2017SS, endyr + 5: objective +1,215). An R0 linkage under a
# hindcast curve does nothing single-species and sets only R_init under predation.

.srr_guard_msg <- function(d, pattern) {
  msg <- tryCatch({ suppressWarnings(suppressMessages(Rceattle:::data_check(d))); "" },
                  error = function(e) conditionMessage(e))
  grepl(pattern, msg)
}

# The bundled data store msmMode = 1 (fit_mod()'s argument overrides it), so set it.
.penalty_data <- function() {
  d <- Rceattle::BS2017SS
  d$msmMode <- 0
  d$srr_fun <- 0L; d$srr_pred_fun <- 2L
  d$srr_hat_styr <- d$styr + 1; d$srr_hat_endyr <- d$endyr
  d
}

test_that("penalty years outside the hindcast are refused", {
  d <- .penalty_data()
  expect_false(.srr_guard_msg(d, "srr_hat_"))

  d1 <- d; d1$srr_hat_endyr <- d$endyr + 1
  expect_true(.srr_guard_msg(d1, "srr_hat_endyr .* must be <= endyr"))

  d2 <- d; d2$srr_hat_styr <- d$styr - 1
  expect_true(.srr_guard_msg(d2, "srr_hat_styr .* must be >= styr"))

  # build_srr(srr_hat_endyr = ) on its own resets srr_pred_fun to 0 and drops the
  # penalty, so the message must point at the model's full build_srr() call.
  expect_true(.srr_guard_msg(d1, "full build_srr\\(\\) call"))
  expect_true(.srr_guard_msg(d2, "full build_srr\\(\\) call"))
  expect_true(.srr_guard_msg(d1, "on its own .*drops the penalty"))
  expect_true(.srr_guard_msg(d2, "on its own .*drops the penalty"))
})

test_that("an empty penalty window, or no penalty at all, is not refused", {
  # A retrospective peel ending before the penalty years has srr_hat_styr > srr_hat_endyr.
  d <- .penalty_data()
  d$srr_hat_styr <- d$endyr; d$srr_hat_endyr <- d$endyr - 5
  expect_false(.srr_guard_msg(d, "srr_hat_"))

  # Without the penalty form the template never reads these years.
  d$srr_pred_fun <- 0L; d$srr_hat_endyr <- d$endyr + 10
  expect_false(.srr_guard_msg(d, "srr_hat_"))
})

test_that("an R0 linkage under a single-species hindcast curve is refused", {
  d <- Rceattle::BS2017SS
  d$msmMode <- 0
  d$srr_fun <- d$srr_pred_fun <- 2L
  d$srr_linkages <- list(R0 = Rceattle::linkage_spec(~ 1))
  expect_true(.srr_guard_msg(d, "R0 linkage has no effect"))

  # A string alias resolves before the comparison ("SingleSpecies" > 0 is TRUE in R).
  d$msmMode <- "SingleSpecies"
  expect_true(.srr_guard_msg(d, "R0 linkage has no effect"))
  d$msmMode <- 0

  # Mean recruitment keeps it, with or without a covariate.
  d$srr_fun <- d$srr_pred_fun <- 0L
  d$srr_linkages <- list(R0 = Rceattle::linkage_spec(~ BTempC))
  expect_false(.srr_guard_msg(d, "R0 linkage|R_init"))
})

test_that("under predation only an intercept-only R0 linkage is accepted", {
  d <- Rceattle::BS2017MS
  d$msmMode <- 1
  d$srr_fun <- d$srr_pred_fun <- 2L
  d$srr_linkages <- list(R0 = Rceattle::linkage_spec(~ 1))
  expect_false(.srr_guard_msg(d, "R0 linkage|R_init"))

  for (f in list(~ BTempC, ~ (1 | Year))) {
    d$srr_linkages <- list(R0 = Rceattle::linkage_spec(f))
    expect_true(.srr_guard_msg(d, "covariate or random effect"), info = deparse(f))
  }

  # A per-species list is checked spec by spec.
  d$srr_linkages <- list(R0 = list(Rceattle::linkage_spec(~ 1, species = 1),
                                   Rceattle::linkage_spec(~ BTempC, species = 2)))
  expect_true(.srr_guard_msg(d, "covariate or random effect"))
})

test_that("under predation an intercept R0 prior lands on R_init", {
  set.seed(123)
  d <- make_msm_test_data()$data_list
  build <- function(inits) suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = inits, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    recFun = build_srr(srr_fun = "BevertonHolt", linkages = list(
      R0 = linkage_spec(~ 1, species = 1,
                        priors = list(`(Intercept)` = prior_lognormal(log(1e4), 0.5))))),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  # Distinct starting levels, so a prior landing on species 2 would not pass.
  p <- build(NULL)$estimated_params
  p$rec_pars[, 1] <- c(8, 10)
  m <- build(p)
  q <- m$quantities
  expect_gt(abs(log(q$R_init[1]) - log(q$R_init[2])), 1)
  expect_true(is.finite(m$obj$fn()))
  expect_equal(sum(q$jnll_comp["Linkage-table priors", ]),
               -dnorm(log(unname(q$R_init[1])), log(1e4) - 0.5^2 / 2, 0.5, log = TRUE),
               tolerance = 1e-8)
})
