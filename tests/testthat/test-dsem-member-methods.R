# The dsem member functions on fit$dsem (issue #220).
#
# build_dsem_objects() returns dsem::dsem(run_model = FALSE), which classes its
# own return "dsem". Left at that, summary() and plot() dispatch to the dsem
# package, which reads the path coefficients out of $obj -- the standalone TMB
# object built BEFORE fitting, so it still holds the starting values. Measured
# on a converged GOApollock DSEM (objective 797.9304), summary(fit$dsem) gave
# AR_BT = 0, BT_to_R = 0, sigmaR1 = 1 against MLEs of 0.0544, -0.2951, 1.3353.
# The fixture below reproduces that: `starts` is what dsem's own method reads,
# `beta_z` what the fit estimated, and they differ on every path.

fake_sem_full <- function() {
  data.frame(
    path      = c("BT -> recdevs1", "recdevs1 <-> recdevs1"),
    lag       = c(1, 0),
    name      = c("BT_to_R", "sigmaR1"),
    start     = c(0, 1),
    parameter = c(1L, 2L),
    first     = c("BT", "recdevs1"),
    second    = c("recdevs1", "recdevs1"),
    direction = c(1, 2),
    stringsAsFactors = FALSE
  )
}

# A fit whose DSEM slot is shaped as build_dsem_objects() leaves it: dsem's
# $obj, $internal and $sem_full, plus the Rceattle subclass. $obj is a stand-in
# that returns the STARTING values, which is all summary.dsem() reads from it.
fake_dsem_fit <- function(beta_z = c(0.5, 0.8), starts = c(0.01, 1)) {
  tsdata <- stats::ts(matrix(0, 3, 2,
                             dimnames = list(NULL, c("recdevs1", "BT"))))
  slot <- structure(
    list(
      sem_full   = fake_sem_full(),
      obj        = list(env = list(parList = function(...) {
        list(beta_z = starts)
      })),
      internal   = list(tsdata = tsdata),
      tmb_inputs = list(data = list(rec_sd_idx = 2L))
    ),
    class = c("Rceattle_dsem", "dsem"))

  fit <- structure(
    list(
      data_list        = list(nspp = 1L, spnames = "spp1", random_rec = TRUE),
      estimated_params = list(dsem_beta_z = beta_z),
      quantities       = list(R_sd = 0.8),
      map              = list(mapList = list(dsem_beta_z = factor(c(1, 2)))),
      sdrep            = NULL,
      dsem             = slot
    ),
    class = "Rceattle")

  # What fit_mod() records once the optimization and sdreport are done.
  Rceattle:::.rce_record_dsem_coefs(fit)
}


testthat::test_that("build_dsem_objects() subclasses its return", {
  testthat::skip_if_not_installed("dsem")
  d <- list(styr = 1990, endyr = 2010, projyr = 2020, nspp = 1L, spnames = "spp1",
            sigma_rec_prior = 1, random_rec = TRUE, proj_mean_rec = FALSE,
            env_data = data.frame(Year = 1990:2020,
                                  BT = as.numeric(scale(seq_len(31)))))
  b <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = "BT -> recdevs1, 1, BT_to_R, 0
                                      recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                               family = "fixed", estimate_projection = FALSE),
    data_list = d)

  testthat::expect_s3_class(b, "Rceattle_dsem")
  # Still a dsem, so anything testing for one keeps working.
  testthat::expect_true(inherits(b, "dsem"))
  testthat::expect_identical(class(b)[1], "Rceattle_dsem")
})


testthat::test_that("summary() reports the estimates, not the starting values", {
  testthat::skip_if_not_installed("dsem")
  fit <- fake_dsem_fit(beta_z = c(0.5, 0.8), starts = c(0.01, 1))

  # The positive case: with the bare dsem class the same object reports its
  # start values, which is the behaviour this subclass exists to displace.
  bare <- fit$dsem
  class(bare) <- "dsem"
  testthat::expect_equal(summary(bare)$Estimate, c(0.01, 1))

  testthat::expect_equal(summary(fit$dsem)$Estimate, c(0.5, 0.8))
})


testthat::test_that("summary() of an unfitted spec says the values are starts", {
  testthat::skip_if_not_installed("dsem")
  b <- fake_dsem_fit()$dsem
  b$path_coefficients <- NULL

  # Not an error: for the object build_dsem_objects() returns, the starting
  # values ARE the answer, and plot.dsem() calls summary() before it draws.
  testthat::expect_message(out <- summary(b), "starting values")
  testthat::expect_equal(out$Estimate, c(0.01, 1))
})


testthat::test_that("plot() labels the path diagram with the estimates", {
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not_installed("igraph")
  fit <- fake_dsem_fit(beta_z = c(0.5, 0.8))

  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())

  # The control, not just the result: dsem's plot.dsem() builds its labels by
  # calling summary(x), so drawing the SAME object under the bare class is what
  # shows the subclass is why the arrows carry MLEs instead of start values.
  bare <- fit$dsem
  class(bare) <- "dsem"
  testthat::expect_equal(igraph::E(plot(bare, edge_label = "value"))$label,
                         c(0.01, 1))
  testthat::expect_equal(igraph::E(plot(fit$dsem, edge_label = "value"))$label,
                         c(0.5, 0.8))
})


testthat::test_that("the methods that cannot be honest refuse, and say what to use", {
  d <- fake_dsem_fit()$dsem

  testthat::expect_error(logLik(d),    "jnll_dsem")
  testthat::expect_error(vcov(d),      "vcov\\(fit\\)")
  testthat::expect_error(residuals(d), "process_residuals")
  testthat::expect_error(predict(d),   "dsem_z_tj")
  testthat::expect_error(simulate(d),  "sim_mod")

  # Each says why, not only what.
  for (f in list(logLik, vcov, residuals, predict, simulate)) {
    testthat::expect_error(f(d), "joint CEATTLE objective")
  }
})


testthat::test_that(".as_rceattle_dsem() prepends once, and adds no 'list'", {
  x <- structure(list(), class = "dsem")
  testthat::expect_identical(class(Rceattle:::.as_rceattle_dsem(x)),
                             c("Rceattle_dsem", "dsem"))
  # Idempotent: a re-passed slot must not stack the class.
  testthat::expect_identical(
    class(Rceattle:::.as_rceattle_dsem(Rceattle:::.as_rceattle_dsem(x))),
    c("Rceattle_dsem", "dsem"))
  # oldClass(), not class(): an unclassed list gains no implicit "list".
  testthat::expect_identical(class(Rceattle:::.as_rceattle_dsem(list())),
                             "Rceattle_dsem")
})


testthat::test_that("re-recording overwrites a table the estimates have moved past", {
  fit <- fake_dsem_fit(beta_z = c(0.5, 0.8))
  testthat::expect_equal(summary(fit$dsem)$Estimate, c(0.5, 0.8))

  # What retrospective() relies on after it swaps a peel's sdreport.
  fit$estimated_params$dsem_beta_z <- c(1.5, 2.5)
  again <- Rceattle:::.rce_record_dsem_coefs(fit)
  testthat::expect_equal(summary(again$dsem)$Estimate, c(1.5, 2.5))
})


testthat::test_that("summary(fit) warns when the slot's table has gone stale", {
  testthat::skip_if_not_installed("dsem")
  fit <- fake_dsem_fit(beta_z = c(0.5, 0.8))

  # A mutation with no re-record -- what model_average() and retrospective()
  # would leave behind if their handling were dropped.
  fit$estimated_params$dsem_beta_z <- c(1.5, 2.5)
  testthat::expect_warning(suppressMessages(summary(fit)), "stale path table")
})


testthat::test_that("model_average() leaves no component's path table behind", {
  testthat::skip_if_not_installed("dsem")
  fit <- fake_dsem_fit(beta_z = c(0.5, 0.8))

  # The shape model_average() produces: one component fit, estimates dropped.
  fit$dsem$path_coefficients <- NULL
  testthat::expect_message(summary(fit$dsem), "starting values")
})


testthat::test_that("print() shows the variables and the table", {
  fit <- fake_dsem_fit()
  out <- utils::capture.output(print(fit$dsem))
  testthat::expect_match(out[1], "DSEM on 2 variables: recdevs1, BT")
  testthat::expect_true(any(grepl("BT_to_R", out)))

  b <- fit$dsem
  b$path_coefficients <- NULL
  testthat::expect_match(paste(utils::capture.output(print(b)), collapse = " "),
                         "model inputs only")
})


# The only code that records the table on a real fit is fit_mod()'s call to
# .rce_record_dsem_coefs(), and it must run AFTER `sdrep` is attached: move it
# above and every Std_Error silently becomes NA while the suite stays green.
testthat::test_that("fit_mod() records the path table, with its standard errors", {
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not_installed("TMB")
  testthat::skip_on_cran()

  data("GOApollock", package = "Rceattle")
  GOApollock$projyr <- 2020
  GOApollock$env_data$ScaledBT <- as.numeric(scale(GOApollock$env_data$BTempC))

  fit <- Rceattle::fit_mod(
    data_list = GOApollock, inits = NULL, file = NULL, estimateMode = 0,
    dsem = build_DSEM(sem = "ScaledBT -> ScaledBT, 1, AR_BT, 0
                             ScaledBT -> recdevs1, 1, BT_to_R, 0
                             recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                      family = "fixed"),
    random_rec = TRUE, msmMode = 0, initMode = 2,
    fit_control = fit_control(phase = TRUE, verbose = 0))

  testthat::expect_s3_class(fit$dsem, "Rceattle_dsem")
  testthat::expect_equal(summary(fit$dsem), summary(fit)$dsem_coefficients)

  # The ordering pin: .rce_record_dsem_coefs() must run after fit_mod()
  # attaches `sdrep`. Move it above and every Std_Error comes back NA. Asserted
  # against the fit's own sdreport so a failure cannot be read as "this fit
  # produced no Hessian".
  testthat::expect_false(is.null(fit$sdrep))
  testthat::expect_false(all(is.na(summary(fit$dsem)$Std_Error)))

  # Not the starting values the slot's $obj still holds, derived from the bare
  # class rather than hard-coded, so a length change cannot pass it vacuously.
  bare <- fit$dsem
  class(bare) <- "dsem"
  starts <- summary(bare)$Estimate
  est <- summary(fit$dsem)$Estimate
  testthat::expect_equal(length(est), length(starts))
  testthat::expect_false(isTRUE(all.equal(est, starts)))
})
