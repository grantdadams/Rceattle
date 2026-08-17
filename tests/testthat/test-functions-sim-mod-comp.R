# Age/length compositions and conditional age-at-length are drawn by the TMB
# model, in RAW bin space (ceattle.cpp, slots 2 and 3).
#
# The load-bearing claim is that drawing raw and letting the refit fold is
# EXACT, not an approximation: tail accumulation merges bins many-to-one before
# the density, and a draw taken there could not be written back at all. It is
# exact because the multinomial and the Dirichlet-multinomial are both closed
# under merging categories. These tests check that against the folded
# distribution rather than assuming it.
testthat::skip_on_cran()

.comp_fixture <- function(nages = 6, N = 200, yng = NULL, old = NULL,
                          dist = NULL, weight = NULL, seed = 123) {
  d <- make_test_data(nyrs = 20, nages = nages, seed = seed)
  d$comp_data$Sample_size <- N
  if (!is.null(yng))    d$fleet_control$Comp_accum_young <- yng
  if (!is.null(old))    d$fleet_control$Comp_accum_old   <- old
  if (!is.null(dist))   d$fleet_control$Comp_distribution <- dist
  if (!is.null(weight)) d$fleet_control$Comp_weights <- weight
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
}

.comp_draws <- function(fit, nrep, seed = 9) {
  cols <- grep("^Comp_", names(fit$data_list$comp_data))
  set.seed(seed)
  replicate(nrep, suppressWarnings(
    as.matrix(Rceattle::sim_mod(fit, simulate = TRUE)$comp_data[, cols])))
}


testthat::test_that("a raw draw folds to the folded multinomial", {
  testthat::skip_if_not_installed("TMB")

  yng <- 2; old <- 4
  fit  <- .comp_fixture(nages = 6, N = 200, yng = yng, old = old)
  cols <- grep("^Comp_", names(fit$data_list$comp_data))
  hat  <- fit$quantities$comp_hat[, seq_along(cols), drop = FALSE]
  N    <- fit$data_list$comp_data$Sample_size

  fold <- function(v) {
    as.numeric(tapply(v, pmin(pmax(seq_along(v), yng), old), sum))
  }

  reps <- .comp_draws(fit, nrep = 800)
  row  <- 1
  raw  <- t(reps[row, , ])
  # Every raw draw is a complete sample: it is only the FOLD that must match.
  testthat::expect_true(all(abs(rowSums(raw) - N[row]) < 1e-8))

  folded <- t(apply(raw, 1, fold))
  p      <- fold(hat[row, ])
  keep   <- which(N[row] * p > 5)
  testthat::expect_gt(length(keep), 1)

  z <- (colMeans(folded) - N[row] * p) /
    (sqrt(N[row] * p * (1 - p)) / sqrt(nrow(folded)))
  testthat::expect_lt(max(abs(z[keep])), 4)
  # And the spread is the folded multinomial's, not the raw one's.
  testthat::expect_equal(apply(folded, 2, stats::sd)[keep],
                         sqrt(N[row] * p * (1 - p))[keep], tolerance = 0.15)
})


testthat::test_that("the Dirichlet-multinomial composition is over-dispersed", {
  testthat::skip_if_not_installed("TMB")

  # theta is exp(Comp_weights) and Comp_weights is ESTIMATED, so a fitted model
  # drives it large, where the Dirichlet-multinomial is a multinomial to machine
  # precision and its expected over-dispersion is 1.0000 -- a check run there
  # passes whatever the simulator does. Hold every parameter at one fit and vary
  # only theta, exactly as the diet harness has to.
  base <- .comp_fixture(dist = "Multinomial", weight = 1)
  hold <- function(dist, w) {
    d <- base$data_list
    d$fleet_control$Comp_distribution <- dist
    inits <- base$estimated_params
    inits$comp_weights[] <- w
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, inits = inits, file = NULL, estimateMode = 3, msmMode = 0,
      random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE))))
  }
  mn <- hold("Multinomial", 1)            # weight 1: draw at the nominal N
  dm <- hold("DirichletMultinomial", 0)   # theta = exp(0) = 1

  cols <- grep("^Comp_", names(mn$data_list$comp_data))
  hat  <- mn$quantities$comp_hat[, seq_along(cols), drop = FALSE]
  N    <- mn$data_list$comp_data$Sample_size[1]
  big  <- which(hat[1, ] > 0.05)
  testthat::expect_gt(length(big), 1)

  sd_mn <- apply(t(.comp_draws(mn, 600)[1, , ])[, big, drop = FALSE], 2, stats::sd)
  sd_dm <- apply(t(.comp_draws(dm, 600)[1, , ])[, big, drop = FALSE], 2, stats::sd)

  expected <- sqrt((N + N * 1) / (1 + N * 1))   # ~1.41 at N = 200, theta = 1
  testthat::expect_gt(expected, 1.3)            # the check has power
  testthat::expect_equal(mean(sd_dm / sd_mn), expected, tolerance = 0.2)
})


testthat::test_that("Comp_weights scales the multinomial draw", {
  testthat::skip_if_not_installed("TMB")

  # The weight multiplies the multinomial log-likelihood, which makes it an
  # effective sample size. Drawing at the nominal N regardless would hand the
  # estimator data of the wrong precision.
  w <- 4
  f1 <- .comp_fixture(weight = 1)
  fw <- .comp_fixture(weight = w)

  cols <- grep("^Comp_", names(f1$data_list$comp_data))
  N    <- f1$data_list$comp_data$Sample_size[1]

  r1 <- t(.comp_draws(f1, 400)[1, , ])
  rw <- t(.comp_draws(fw, 400)[1, , ])
  testthat::expect_equal(unname(rowSums(r1)[1]), N)
  testthat::expect_equal(unname(rowSums(rw)[1]), N * w)

  hat <- f1$quantities$comp_hat[1, seq_along(cols)]
  big <- which(hat > 0.05)
  # Proportions, so the extra sample size shows as a sqrt(w) tighter spread.
  testthat::expect_equal(
    mean(apply(rw[, big, drop = FALSE] / (N * w), 2, stats::sd) /
         apply(r1[, big, drop = FALSE] / N, 2, stats::sd)),
    1 / sqrt(w), tolerance = 0.15)
})


testthat::test_that("an empty composition row is reported, not passed through", {
  testthat::skip_if_not_installed("TMB")

  # A row sums to zero when the model predicts no composition for it -- a fleet
  # with no catch that year, or one switched off -- or when Sample_size times the
  # composition weight rounds away. run_mse() reads a zero row as "not sampled"
  # and drops the sample size with it, so the row must come back zero rather than
  # holding its original proportions (which sum to 1 and look entirely
  # plausible). Say so rather than let a self_test() score a replicate whose
  # data quietly went missing.
  testthat::expect_warning(
    .sim_warn_empty_comp(matrix(c(0, 0, 1, 2), nrow = 2, byrow = TRUE),
                         c(1, 2), "composition"),
    "came back empty")
  testthat::expect_no_warning(
    .sim_warn_empty_comp(matrix(c(3, 1, 1, 2), nrow = 2, byrow = TRUE),
                         c(1, 2), "composition"))
})


testthat::test_that("CAAL rows the model predicts nothing for come back empty", {
  testthat::skip_if_not_installed("TMB")

  # NOTE what this does and does not cover. No fixture in this repo produces a
  # non-zero caal_hat: make_test_data()'s CAAL rows carry Sample_size 1 with no
  # predicted composition, make_msm_test_data()'s carry Sample_size 0, and none
  # of the bundled data sets has CAAL at all. So the CAAL *draw* is still
  # unexercised, and saying so is more useful than a test that looks like
  # coverage. What is checked here is the write-back contract that the draw
  # shares with the marginal composition: a row with no prediction comes back
  # zero rather than holding its own values, which is what run_mse() reads as
  # "not sampled".
  d <- make_test_data(nyrs = 20, nages = 6, seed = 123, growth = "vonBertalanffy")
  testthat::expect_gt(nrow(d$caal_data), 0)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  testthat::expect_equal(sum(rowSums(fit$quantities$caal_hat) > 0), 0)

  cols <- grep("^CAAL_", names(fit$data_list$caal_data))
  set.seed(6)
  sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))
  got <- as.matrix(sim$caal_data[, cols])

  testthat::expect_equal(nrow(got), nrow(fit$data_list$caal_data))
  testthat::expect_false(anyNA(got))
  testthat::expect_true(all(got >= 0))
  # Unpredicted rows are zeroed, not left holding their original proportions.
  testthat::expect_true(all(rowSums(got) == 0))
})


testthat::test_that("simulate = FALSE returns the predicted composition and draws nothing", {
  testthat::skip_if_not_installed("TMB")

  fit  <- .comp_fixture()
  cols <- grep("^Comp_", names(fit$data_list$comp_data))

  set.seed(2); before <- .Random.seed
  sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = FALSE))
  testthat::expect_identical(.Random.seed, before)
  testthat::expect_equal(unname(as.matrix(sim$comp_data[, cols])),
                         unname(fit$quantities$comp_hat[, seq_along(cols)]),
                         tolerance = 1e-10)
})


testthat::test_that("a joint-sex draw folds correctly within each sex block", {
  testthat::skip_if_not_installed("TMB")

  # The structural case. For a joint-sex composition (Sex = 3) the row holds two
  # sex blocks end to end, and the C++ folds the tails WITHIN each block so a
  # fold never crosses the sex boundary -- while the draw is taken across the
  # whole 2*nbins vector at once. That is still exact, because merging categories
  # is merging categories wherever the boundaries fall, but it is the case most
  # likely to be got wrong, so pin it rather than reason about it.
  utils::data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- GOA2018SS
  # Fleet 9 carries joint-sex AGE comps for species 2 (21 bins). Accumulation
  # bins are ordinals on the fleet's OWN dimension, so setting them globally
  # would exceed the bin count of every narrower fleet and data_check() would --
  # correctly -- reject the model. Set them on this fleet alone.
  target_flt <- 9L
  yng <- 4; old <- 15
  d$fleet_control$Comp_accum_young[d$fleet_control$Fleet_code == target_flt] <- yng
  d$fleet_control$Comp_accum_old[d$fleet_control$Fleet_code == target_flt]   <- old
  # This fleet's fitted composition weight is well below 1, so at its own sample
  # size the draw places about one observation per folded bin -- too thin for a
  # per-bin check. Raise the sample size rather than the tolerance.
  d$comp_data$Sample_size[d$comp_data$Fleet_code == target_flt] <- 2000

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))

  cd  <- fit$data_list$comp_data
  row <- which(cd$Sex == 3 & cd$Fleet_code == target_flt &
                 cd$Year > 0 & cd$Year <= fit$data_list$endyr &
                 cd$Sample_size > 0)[1]
  testthat::skip_if(is.na(row), "no fitted joint-sex composition row")

  sp <- cd$Species[row]
  nb <- if (cd$Age0_Length1[row] == 0) fit$data_list$nages[sp] else
    fit$data_list$nlengths[sp]
  testthat::expect_gte(nb, old)
  cols  <- grep("^Comp_", names(cd))
  hat   <- fit$quantities$comp_hat[row, seq_along(cols)]
  N     <- cd$Sample_size[row]

  # Fold each sex block on its own, exactly as ceattle.cpp does.
  fold_blocks <- function(v) {
    unlist(lapply(0:1, function(b) {
      seg <- v[b * nb + seq_len(nb)]
      as.numeric(tapply(seg, pmin(pmax(seq_len(nb), yng), old), sum))
    }))
  }

  set.seed(21)
  reps <- replicate(600, suppressWarnings(
    as.numeric(Rceattle::sim_mod(fit, simulate = TRUE)$comp_data[row, cols])))
  raw <- t(reps)   # replicates x bins

  # The draw is taken at Sample_size * Comp_weights, and this fleet's weight is
  # fitted rather than 1, so use the total actually drawn -- the fold claim is
  # about the distribution GIVEN that total. Every replicate must place the same
  # number, which is what makes the comparison well posed.
  Nd <- rowSums(raw)[1]
  testthat::expect_true(all(abs(rowSums(raw) - Nd) < 1e-8))
  testthat::expect_gt(Nd, 0)

  # Each sex block gets its own share; the split must not drift.
  blk1 <- rowSums(raw[, seq_len(nb), drop = FALSE])
  testthat::expect_equal(mean(blk1) / Nd, sum(hat[seq_len(nb)]), tolerance = 0.05)

  folded <- t(apply(raw, 1, fold_blocks))
  p      <- fold_blocks(hat)
  p      <- p / sum(p)
  keep   <- which(Nd * p > 5)
  testthat::expect_gt(length(keep), 2)

  z <- (colMeans(folded) - Nd * p) / (sqrt(Nd * p * (1 - p)) / sqrt(nrow(folded)))
  testthat::expect_lt(max(abs(z[keep])), 4)
  testthat::expect_equal(apply(folded, 2, stats::sd)[keep],
                         sqrt(Nd * p * (1 - p))[keep], tolerance = 0.2)
})
