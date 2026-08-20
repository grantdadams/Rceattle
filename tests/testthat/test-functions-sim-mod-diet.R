# Stomach contents are drawn by the TMB model's SIMULATE block, beside the diet
# density (ceattle.cpp, section 13.2). Before that they were the one observation
# type never simulated, so a multispecies self_test() refit against the same
# stomachs every replicate and recovered suitability from data that never varied.
#
# The check that matters is the ROUND TRIP: take the simulated stomach
# proportions, put them through the transform a refit applies (append the "other
# prey" balance, add comp_offset, renormalize), and compare their mean to the
# proportions the density predicts. That is what the estimator actually sees.
#
# Comparing the RAW simulated proportions to the offset-normalized prediction
# instead -- the obvious-looking check -- cannot detect a draw that is centred on
# the wrong one of the two, because at the default comp_offset the two targets
# differ by only ~0.02%. Raise comp_offset to give the test power.
testthat::skip_on_cran()

# Per stomach, in diet_data row order: the prey proportions the density predicts
# (offset and renormalized), and the same transform applied to a simulated row.
.diet_transform <- function(mod, values, offset = 1e-5) {
  re  <- Rceattle::rearrange_data(mod$data_list)
  sid <- re$stomach_id
  out <- rep(NA_real_, length(values))
  for (i in unique(sid)) {
    idx <- which(sid == i)
    p <- c(values[idx], max(1 - sum(values[idx]), 0))
    p <- (p + offset) / sum(p + offset)
    out[idx] <- p[seq_along(idx)]
  }
  out
}

# What the density predicts for each prey bin.
.diet_pred <- function(mod, offset = 1e-5) {
  .diet_transform(mod, mod$quantities$diet_hat[, 2], offset = offset)
}

# What the refit would see, given a simulated table.
.diet_seen <- function(mod, sim, offset = 1e-5) {
  .diet_transform(mod, sim, offset = offset)
}

# Fit BS2017MS once with an estimated suitability, so the diet likelihood is
# active and diet_hat is concentrated enough to measure a spread on. Then hold
# every parameter there and vary only the composition family and its
# concentration, so those are the single thing that differs between cases.
# NOTE the parameter `diet_comp_weights` means different things by family: for
# the multinomial it is used raw, as a likelihood weight / effective-sample-size
# multiplier (1 = neutral), and for the Dirichlet-multinomial it is log theta
# (exp() is taken). Passing log(1) = 0 for a multinomial sets the weight to
# ZERO, which drops the data from the likelihood and from the draw.
.diet_sim_env <- new.env(parent = emptyenv())
.diet_sim_fit <- function(family, weight_par) {
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
  inits$diet_comp_weights[] <- weight_par
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

  mod <- .diet_sim_fit("Multinomial", 1)
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


testthat::test_that("what a refit sees is centred on what the density predicts", {
  testthat::skip_if_not_installed("TMB")

  # Run at a raised comp_offset as well as the default. Drawing from the
  # offset-and-renormalized prediction rather than the raw one applies the offset
  # twice over the round trip, which is 0.02% at the default and over 20% at
  # 1e-2 -- invisible at the default, unmissable here.
  for (fam in c("Multinomial", "DirichletMultinomial")) {
    mod <- .diet_sim_fit(fam, if (fam == "Multinomial") 1 else log(1))
    for (off in c(1e-5, 1e-2)) {
      pred <- .diet_pred(mod, offset = off)
      reps <- .diet_reps(mod, nrep = 150)
      seen <- apply(reps, 2, function(x) .diet_seen(mod, x, offset = off))
      keep <- which(pred > 1e-4)
      testthat::expect_gt(length(keep), 50)
      # 0.04 sits well above the Monte Carlo noise of this aggregate and well
      # below the +23% the pre-fix draw carried at comp_offset = 1e-2.
      testthat::expect_equal(sum(rowMeans(seen)[keep]) / sum(pred[keep]), 1,
                             tolerance = 0.04,
                             label = paste(fam, "comp_offset", off))
    }
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
  mn <- .diet_sim_fit("Multinomial", 1)
  dm <- .diet_sim_fit("DirichletMultinomial", log(0.1))

  tgt <- .diet_pred(mn)
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


testthat::test_that("the multinomial diet draw honours Diet_comp_weights", {
  testthat::skip_if_not_installed("TMB")

  # For diet_ll_type == 0 the weight multiplies the log-likelihood, which makes
  # it an effective sample size. Drawing at the nominal N regardless would give
  # a down-weighted predator stomachs the likelihood treats as w times less
  # informative than they are, and self_test() would report recovery against
  # data of the wrong precision. (Note the weight lives in the PARAMETER
  # diet_comp_weights, not the data column, which fit_mod() ignores when inits
  # are supplied.)
  w <- 4
  mn1 <- .diet_sim_fit("Multinomial", 1)
  mnw <- .diet_sim_fit("Multinomial", w)
  testthat::expect_equal(as.numeric(mnw$estimated_params$diet_comp_weights)[1], w)

  tgt <- .diet_pred(mn1)
  big <- which(tgt > 0.02)
  testthat::expect_gt(length(big), 20)

  sd1 <- apply(.diet_reps(mn1, nrep = 150)[big, , drop = FALSE], 1, stats::sd)
  sdw <- apply(.diet_reps(mnw, nrep = 150)[big, , drop = FALSE], 1, stats::sd)

  # Weighting by w shrinks the spread by sqrt(w).
  testthat::expect_equal(mean(sdw / sd1), 1 / sqrt(w), tolerance = 0.12)
})


testthat::test_that("a diet row whose weighted sample size rounds to zero comes back empty", {
  testthat::skip_if_not_installed("TMB")

  # Regression: the write-back used to be gated on the draw having placed
  # something, so a stomach whose N_s * Diet_comp_weights fell below one
  # observation silently kept its OBSERVED proportions. A self_test() then
  # refitted against the real diet on those rows and read suitability recovery
  # as better than it was. The per-predator "frozen" warning could not see it --
  # a predator with some stomachs drawn and some rounded away is not frozen.
  # Compositions have always zeroed such a row; diet now does too.
  mod <- .diet_sim_fit("Multinomial", 1e-4)

  obs <- mod$data_list$diet_data$Stomach_proportion_by_weight
  sim <- suppressWarnings(
    Rceattle::sim_mod(mod, simulate = TRUE))$diet_data$Stomach_proportion_by_weight

  # Every row the model fits is emptied, not carried over.
  had <- obs > 0
  testthat::expect_gt(sum(had), 0)
  testthat::expect_true(all(sim[had] == 0))

  # And it is announced rather than silent.
  testthat::expect_warning(Rceattle::sim_mod(mod, simulate = TRUE),
                           "came back empty")
})
