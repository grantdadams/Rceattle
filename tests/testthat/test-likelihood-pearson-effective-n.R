# Pearson residuals must divide by the effective sample size the composition
# likelihood actually used, not by the raw input sample size.
#
# Two things make those differ. Under a multinomial (and the AFSC
# pseudo-likelihood) the weight column multiplies the log-likelihood, so it is an
# effective sample size -- ceattle.cpp:3758 draws the simulated composition at
# `n_nom * comp_weights(flt)`, which states the assumed variance outright. Under
# a Dirichlet-multinomial the proportions are overdispersed by
# (n + conc)/(1 + conc). Before v5.29.0 both were ignored: every composition
# residual divided by the multinomial variance at the raw Sample_size, so a
# weighted or DM fleet's residuals were inflated and read as systematic misfit.
#
# These are Monte Carlo self-consistency checks, deliberately not restatements of
# the formula: data are drawn under the family the likelihood assumes and the
# residuals must come back standard normal. A test that recomputed the closed
# form in R and compared it to the same closed form would pass on a wrong
# variance.

set.seed(4321)

# Draw one Dirichlet vector.
rdirich <- function(a) { g <- stats::rgamma(length(a), a, 1); g / sum(g) }

# Per-bin sd of the Pearson residual across replicates. Under a correct variance
# this is 1: the residual is standardized marginally, bin by bin.
resid_sd <- function(draw, K, nrep) {
  r <- vapply(seq_len(nrep), function(i) draw(), numeric(K))
  apply(r, 1L, stats::sd)
}


test_that("a weighted multinomial fleet's Pearson residuals are standard normal", {
  K <- 8L; N <- 150; w <- 3.5
  p <- c(0.05, 0.10, 0.18, 0.22, 0.20, 0.13, 0.08, 0.04)
  row <- rep(1L, K)

  # The likelihood behaves as a multinomial at w * N, so draw there.
  n_eff <- round(w * N)
  s <- resid_sd(function() {
    obs <- stats::rmultinom(1, n_eff, p)[, 1] / n_eff
    .rce_comp_pearson(obs, p, N, row, family = rep(0L, K), weight = rep(w, K),
                      offset = 0)
  }, K, 4000)

  expect_equal(unname(s), rep(1, K), tolerance = 0.06)
})


test_that("ignoring the weight rescales a multinomial residual by 1/sqrt(w)", {
  # The pre-5.29.0 behaviour, kept as the regression anchor. Dividing by the raw
  # N assumes a variance w times too large, so the old residual is 1/sqrt(w)
  # times the right one: too SMALL on an upweighted fleet (w > 1), too large on
  # a downweighted one. reweight_comps() commonly tunes below 1, so both
  # directions occur in practice.
  K <- 8L; N <- 150
  p <- c(0.05, 0.10, 0.18, 0.22, 0.20, 0.13, 0.08, 0.04)

  for (w in c(3.5, 0.4)) {
    n_eff <- round(w * N)
    obs   <- stats::rmultinom(1, n_eff, p)[, 1] / n_eff
    correct <- as.numeric(
      .rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L, K), rep(w, K), 0))
    old <- (obs - p) / sqrt(p * (1 - p) / N)
    # Compared as a product rather than a ratio: a draw can land exactly on the
    # fitted proportion, and 0/0 is not informative about the scaling.
    expect_equal(old, correct / sqrt(w), tolerance = 1e-8,
                 info = paste("w =", w))
  }
})


test_that("the AFSC pseudo-likelihood is treated as a multinomial", {
  K <- 6L; N <- 90; w <- 2
  p <- c(0.1, 0.2, 0.3, 0.2, 0.15, 0.05)
  obs <- stats::rmultinom(1, round(w * N), p)[, 1] / round(w * N)
  expect_equal(
    .rce_comp_pearson(obs, p, N, rep(1L, K), rep(-1L, K), rep(w, K), 0),
    .rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L,  K), rep(w, K), 0))
})


test_that("a Dirichlet-multinomial fleet's Pearson residuals are standard normal", {
  K <- 10L; N <- 200
  p <- c(0.04, 0.08, 0.14, 0.18, 0.17, 0.14, 0.10, 0.07, 0.05, 0.03)
  row <- rep(1L, K)

  for (theta in c(0.05, 0.5)) {
    lw <- rep(log(theta), K)     # the column is a LOG concentration under a DM
    s <- resid_sd(function() {
      # alpha = N * p * theta, then counts at N -- the cpp's construction.
      obs <- stats::rmultinom(1, N, rdirich(N * p * theta))[, 1] / N
      .rce_comp_pearson(obs, p, N, row, family = rep(1L, K), weight = lw,
                        offset = 0)
    }, K, 4000)
    expect_equal(unname(s), rep(1, K), tolerance = 0.07,
                 info = paste("theta =", theta))
  }
})


test_that("a large Dirichlet-multinomial concentration recovers the multinomial", {
  K <- 6L; N <- 120
  p <- c(0.1, 0.2, 0.3, 0.2, 0.15, 0.05)
  obs <- stats::rmultinom(1, N, p)[, 1] / N

  dm <- .rce_comp_pearson(obs, p, N, rep(1L, K), rep(1L, K), rep(log(1e8), K), 0)
  mn <- .rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L, K), rep(1, K), 0)
  expect_equal(dm, mn, tolerance = 1e-6)
})


test_that("an unweighted multinomial with no offset reproduces the old residual", {
  # Weight 1, offset 0 is the one configuration where the pre-5.29.0 formula was
  # right; it must stay bit-identical so the change is attributable.
  K <- 7L; N <- 80
  p   <- c(0.05, 0.15, 0.25, 0.25, 0.15, 0.10, 0.05)
  obs <- stats::rmultinom(1, N, p)[, 1] / N
  expect_identical(
    as.numeric(.rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L, K),
                                 rep(1, K), 0)),
    (obs - p) / sqrt(p * (1 - p) / N))
})


test_that("the composition offset is carried onto the likelihood's own scale", {
  # comp_offset (default 1e-5) is added to both proportions before the density,
  # so the fitted proportions no longer sum to one. The residual is taken on the
  # renormalized scale the density sees, which matters most in the smallest bins.
  K <- 10L; N <- 100; off <- 1e-5
  p   <- c(0.5, 0.4, rep(0.1 / 8, 8))
  obs <- stats::rmultinom(1, N, p)[, 1] / N

  S     <- sum(p + off)
  want  <- ((obs + off) / S - (p + off) / S) /
    sqrt(((p + off) / S) * (1 - (p + off) / S) / (N * S))
  expect_equal(
    as.numeric(.rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L, K),
                                 rep(1, K), off)),
    want, tolerance = 1e-12)

  # And it is not a no-op: the smallest bins move by ~0.1% here.
  no_off <- .rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L, K), rep(1, K), 0)
  with_off <- .rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L, K), rep(1, K), off)
  expect_false(isTRUE(all.equal(no_off, with_off, tolerance = 1e-6)))
})


test_that("each observation's total is summed over its own bins", {
  # Two observations in one call must not share a total; a fleet with ragged or
  # tail-folded rows depends on this.
  K <- 5L
  p1 <- c(0.2, 0.2, 0.2, 0.2, 0.2); p2 <- c(0.4, 0.3, 0.2, 0.05, 0.05)
  o1 <- c(0.3, 0.1, 0.2, 0.2, 0.2);  o2 <- c(0.5, 0.2, 0.2, 0.05, 0.05)

  both <- as.numeric(.rce_comp_pearson(c(o1, o2), c(p1, p2),
                                       rep(c(60, 90), each = K),
                                       rep(1:2, each = K), rep(0L, 2 * K),
                                       rep(c(2, 3), each = K), 1e-5))
  sep1 <- as.numeric(.rce_comp_pearson(o1, p1, 60, rep(1L, K), rep(0L, K),
                                       rep(2, K), 1e-5))
  sep2 <- as.numeric(.rce_comp_pearson(o2, p2, 90, rep(1L, K), rep(0L, K),
                                       rep(3, K), 1e-5))
  expect_equal(both, c(sep1, sep2))
})
