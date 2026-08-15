# Stomach contents are drawn by the TMB model's SIMULATE block, beside the diet
# density (ceattle.cpp, section 13.2). Before that they were the one observation
# type never simulated, so a multispecies self_test() refit against the same
# stomachs every replicate and recovered suitability from data that never varied.
#
# The draw is centred on the proportions the DENSITY uses, which are not
# diet_hat: the likelihood appends an "other prey" bin, adds comp_offset and
# renormalizes first. Checking against diet_hat instead looks like a large bias
# on rows where diet_hat is near the offset.
testthat::skip_on_cran()

# The proportions the density is taken on, per stomach, in diet_data row order.
.diet_target <- function(mod) {
  re  <- Rceattle::rearrange_data(mod$data_list)
  sid <- re$stomach_id
  hat <- mod$quantities$diet_hat[, 2]
  out <- rep(NA_real_, length(hat))
  for (i in unique(sid)) {
    idx <- which(sid == i)
    p <- c(hat[idx], 1 - sum(hat[idx]))
    p <- (p + 1e-5) / sum(p + 1e-5)
    out[idx] <- p[seq_along(idx)]
  }
  out
}

# Fit BS2017MS once with an estimated suitability, so the diet likelihood is
# active and diet_hat is concentrated enough to measure a spread on. Then hold
# every parameter there and vary only the composition family and its
# concentration, so those are the single thing that differs between cases.
.diet_sim_env <- new.env(parent = emptyenv())
.diet_sim_fit <- function(family, log_theta) {
  if (is.null(.diet_sim_env$base)) {
    utils::data("BS2017SS", package = "Rceattle", envir = environment())
    utils::data("BS2017MS", package = "Rceattle", envir = environment())
    ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      BS2017SS, inits = NULL, msmMode = 0, estimateMode = 1,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
    .diet_sim_env$dat <- BS2017MS
    .diet_sim_env$base <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      BS2017MS, inits = ss$estimated_params, msmMode = 1, suitMode = 4,
      estimateMode = 1, niter = 3,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  }
  d <- .diet_sim_env$dat
  d$Diet_distribution <- rep(family, d$nspp)
  inits <- .diet_sim_env$base$estimated_params
  inits$diet_comp_weights[] <- log_theta
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, inits = inits, msmMode = 1, suitMode = 4, estimateMode = 3, niter = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
}

.diet_reps <- function(mod, nrep, seed = 7) {
  set.seed(seed)
  replicate(nrep, suppressWarnings(
    Rceattle::sim_mod(mod, simulate = TRUE)$diet_data$Stomach_proportion_by_weight))
}


testthat::test_that("simulated stomach contents are a valid diet composition", {
  testthat::skip_if_not_installed("TMB")

  mod <- .diet_sim_fit("Multinomial", log(1))
  set.seed(2)
  sim <- suppressWarnings(Rceattle::sim_mod(mod, simulate = TRUE))
  x <- sim$diet_data$Stomach_proportion_by_weight

  testthat::expect_equal(nrow(sim$diet_data), nrow(mod$data_list$diet_data))
  testthat::expect_false(anyNA(x))
  testthat::expect_true(all(x >= 0))

  # Prey proportions are a share of the stomach; the balance is "other prey",
  # which is not stored and is recomputed on the next fit. Tested at exactly 1,
  # matching data_check(): a looser bound here would pass tables the pipeline
  # rejects, and a rejected refit surfaces as a failed self_test() replicate
  # rather than as bad data.
  sid <- Rceattle::rearrange_data(mod$data_list)$stomach_id
  testthat::expect_true(all(tapply(x, sid, sum) <= 1))

  # The draw actually varies -- diet used to pass through untouched.
  testthat::expect_gt(sum(abs(x - mod$data_list$diet_data$Stomach_proportion_by_weight) > 1e-12), 0)
})


testthat::test_that("a model without an estimated suitability leaves diet alone", {
  testthat::skip_if_not_installed("TMB")

  # suitMode = 0 is empirical suitability: the diet likelihood is not evaluated,
  # nothing is predicted, and so nothing should be redrawn.
  utils::data("BS2017SS", package = "Rceattle", envir = environment())
  ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    BS2017SS, inits = NULL, msmMode = 0, estimateMode = 1,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  set.seed(3)
  sim <- suppressWarnings(Rceattle::sim_mod(ss, simulate = TRUE))
  testthat::expect_identical(sim$diet_data, ss$data_list$diet_data)
})


testthat::test_that("the diet draw is centred where the density expects", {
  testthat::skip_if_not_installed("TMB")

  for (fam in c("Multinomial", "DirichletMultinomial")) {
    mod  <- .diet_sim_fit(fam, log(1))
    tgt  <- .diet_target(mod)
    reps <- .diet_reps(mod, nrep = 150)
    keep <- which(tgt > 1e-4)
    testthat::expect_gt(length(keep), 50)
    # Aggregate rather than a mean of per-row ratios: a ratio on a bin whose
    # target is near the offset is dominated by Monte Carlo noise.
    testthat::expect_equal(sum(rowMeans(reps)[keep]) / sum(tgt[keep]), 1,
                           tolerance = 0.02)
  }
})


testthat::test_that("the Dirichlet-multinomial diet draw is over-dispersed", {
  testthat::skip_if_not_installed("TMB")

  # A Dirichlet-multinomial at theta = 0.1 has sd sqrt((N + N*theta)/(1 + N*theta))
  # times the multinomial's -- about 2.7 at this sample size. Comparing the two
  # draws against each other needs no theoretical constant, and it is the check
  # that has power: at the concentration this model fits to (theta ~ 1e7) the
  # Dirichlet-multinomial is numerically a multinomial, so a test run there would
  # pass whatever the simulator did.
  mn <- .diet_sim_fit("Multinomial", log(1))
  dm <- .diet_sim_fit("DirichletMultinomial", log(0.1))

  tgt <- .diet_target(mn)
  big <- which(tgt > 0.05)
  testthat::expect_gt(length(big), 20)

  sd_mn <- apply(.diet_reps(mn, nrep = 150)[big, , drop = FALSE], 1, stats::sd)
  sd_dm <- apply(.diet_reps(dm, nrep = 150)[big, , drop = FALSE], 1, stats::sd)

  N <- as.numeric(Rceattle::rearrange_data(mn$data_list)$diet_obs[1, 1])
  expected <- sqrt((N + N * 0.1) / (1 + N * 0.1))
  testthat::expect_equal(mean(sd_dm / sd_mn), expected, tolerance = 0.15)
})


testthat::test_that("a fractional Sample_size still yields a valid composition", {
  testthat::skip_if_not_installed("TMB")

  # The draw places a whole number of observations, so dividing the counts by a
  # fractional Sample_size would give proportions summing to round(N)/N > 1 --
  # which data_check() rejects, and self_test() then reports as a replicate that
  # failed to converge. Divide by the number actually placed instead.
  utils::data("BS2017SS", package = "Rceattle", envir = environment())
  utils::data("BS2017MS", package = "Rceattle", envir = environment())
  ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    BS2017SS, inits = NULL, msmMode = 0, estimateMode = 1,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  d <- BS2017MS
  d$diet_data$Sample_size <- 20.6
  mod <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, inits = ss$estimated_params, msmMode = 1, suitMode = 4,
    estimateMode = 1, niter = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  sid <- Rceattle::rearrange_data(mod$data_list)$stomach_id
  set.seed(12)
  worst <- 0
  for (i in 1:15) {
    x <- suppressWarnings(Rceattle::sim_mod(mod, simulate = TRUE))$diet_data$Stomach_proportion_by_weight
    worst <- max(worst, max(tapply(x, sid, sum)))
  }
  testthat::expect_lte(worst, 1)
})
