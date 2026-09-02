# The fold claim: a draw taken in RAW bin space folds to the folded
# distribution EXACTLY, not approximately. Tail accumulation merges bins
# many-to-one before the density, and a draw taken there could not be written
# back at all. It is exact because the multinomial and the
# Dirichlet-multinomial are both closed under merging categories. These tests
# check that against the folded distribution rather than assuming it.
#
# Split out of test-functions-sim-mod-comp.R, which was 27% of the coverage
# suite in one file. Shared fixtures: helpers-comp-sim.R.
testthat::skip_on_cran()

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
