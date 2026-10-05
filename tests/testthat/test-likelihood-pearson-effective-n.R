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


test_that("a multinomial fleet switched off has no Pearson residual", {
  # Comp_weights = 0 drops a multinomial fleet from the likelihood, so it has no
  # effective sample size. Dividing anyway gave sd = Inf and a residual of 0 in
  # every bin, which plot_comp() drew as a perfect fit.
  K <- 6L; N <- 90
  p   <- c(0.1, 0.2, 0.3, 0.2, 0.15, 0.05)
  obs <- stats::rmultinom(1, N, p)[, 1] / N

  r <- .rce_comp_pearson(obs, p, N, rep(1L, K), rep(0L, K), rep(0, K), 0)
  expect_true(all(is.na(as.numeric(r))))
  expect_true(all(is.na(attr(r, "sd"))))

  # A Dirichlet-multinomial reads that column as a log, so 0 is a weight of 1.
  dm <- .rce_comp_pearson(obs, p, N, rep(1L, K), rep(1L, K), rep(0, K), 0)
  expect_true(all(is.finite(as.numeric(dm))))
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


test_that("a Dirichlet-multinomial diet residual is standard normal", {
  # Diet is the third construction: the stomach's proportions are renormalized
  # over prey plus an "other prey" balance BEFORE being scaled by the stomach
  # sample size, so the count total is N_s and the concentration total is
  # N_s * theta -- no second factor, unlike comp and CAAL. Drawn here on that
  # scale, so a residual standardized by any other variance fails.
  K <- 5L; N <- 60; theta <- 0.6
  p <- c(0.30, 0.22, 0.18, 0.12, 0.08)   # prey bins; 0.10 goes to "other prey"

  s <- resid_sd(function() {
    obs <- stats::rmultinom(1, N, rdirich(N * c(p, 0.10) * theta))[, 1] / N
    p_o <- obs[seq_len(K)]
    p_h <- p
    # The residual frame carries prey rows only; the balance is rebuilt from them.
    S_o <- sum(p_o) + (1 - min(sum(p_o), 1))
    S_h <- sum(p_h) + max(1 - sum(p_h), 1e-5)
    .pearson_proportion(p_o / S_o, p_h / S_h, N, conc = N * theta)
  }, K, 4000)

  expect_equal(unname(s), rep(1, K), tolerance = 0.07)
})


test_that("a Dirichlet-multinomial with no recoverable weight falls back to the multinomial", {
  # The fallback has to BE the multinomial, not the Dirichlet-multinomial at a
  # substituted weight. `weight` is a log under a DM, so a stand-in of 1 would
  # mean a concentration of e -- an overdispersion factor near 1.4 at n = 1000,
  # which looks like a considered choice and is not one. The warning says
  # "multinomial variance", so the code must deliver that.
  K <- 5L; N <- 1000
  p <- c(0.30, 0.25, 0.20, 0.15, 0.10)
  o <- c(0.32, 0.24, 0.19, 0.16, 0.09)

  expect_warning(
    .rce_resolve_family_weight(rep(1L, K), rep(NA_real_, K), 1:K, "fleet"),
    "multinomial variance")
  fw <- suppressWarnings(
    .rce_resolve_family_weight(rep(1L, K), rep(NA_real_, K), 1:K, "fleet"))
  expect_identical(fw$family, rep(0L, K))

  got  <- .rce_comp_pearson(o, p, N, rep(1L, K), fw$family, fw$weight, 0)
  want <- .rce_comp_pearson(o, p, N, rep(1L, K), rep(0L, K), rep(1, K), 0)
  expect_equal(as.numeric(got), as.numeric(want))

  # And it is materially different from the theta = e reading it replaces.
  dm_e <- .rce_comp_pearson(o, p, N, rep(1L, K), rep(1L, K), rep(1, K), 0)
  expect_gt(max(abs(as.numeric(dm_e) / as.numeric(want) - 1)), 0.1)

  # A row whose weight IS recoverable is untouched by the fallback.
  ok <- .rce_resolve_family_weight(c(1L, 0L), c(log(2), 3), 1:2, "fleet")
  expect_identical(ok$family, c(1L, 0L))
  expect_identical(ok$weight, c(log(2), 3))
})


test_that("truncating residuals onto the bubble scale keeps sign and reports once", {
  # scale_size_continuous() sets an out-of-bounds size to NA and drops the
  # point, so the largest residuals would be the ones to vanish. The warning
  # names the count and the maximum, not every value: dividing by the effective
  # sample size scales residuals by sqrt(w), which can put thousands past the
  # cap on one panel.
  x   <- c(-9, -6, -1.5, 0, 2, 6, 7.25, NA, NaN)
  got <- suppressWarnings(.rce_truncate_resid(x, "Panel"))
  expect_warning(.rce_truncate_resid(x, "Panel"), "2 residual\\(s\\) beyond")
  expect_warning(.rce_truncate_resid(x, "Panel"), "largest \\|residual\\| 9.00")

  expect_equal(got[1:7], c(-6, -6, -1.5, 0, 2, 6, 6))
  # NA and NaN are not residuals to truncate and must survive as themselves.
  expect_true(is.na(got[8]) && !is.nan(got[8]))
  expect_true(is.nan(got[9]))

  # An infinite residual is still a residual and must land on the cap, not be
  # dropped by the scale.
  expect_equal(suppressWarnings(.rce_truncate_resid(c(Inf, -Inf), "Panel")),
               c(6, -6))

  # Nothing past the cap is silent.
  expect_silent(.rce_truncate_resid(c(-6, 0, 5.9, 6), "Panel"))
})


test_that("a fit saved before a switch rename still resolves its family", {
  # switch_check() upgrades a deprecated spelling at build time, but a fit SAVED
  # before the rename carries the old name, and residuals() runs on the saved
  # object. Finding nothing there would fall back to the schema default -- the
  # multinomial -- which is the same silent wrong-family failure the canonical
  # lookup exists to prevent.
  fc <- data.frame(Fleet_code = 1L, Comp_loglike = 1L, Comp_weights = log(3))
  obj <- structure(list(data_list = list(fleet_control = fc),
                        estimated_params = list()), class = "Rceattle")
  expect_identical(
    .rce_comp_family(obj, 1L, "Comp_distribution", "Comp_weights",
                     "comp_weights")$family, 1L)

  # The canonical spelling wins when both are somehow present.
  fc$Comp_distribution <- 0L
  obj$data_list$fleet_control <- fc
  expect_identical(
    .rce_comp_family(obj, 1L, "Comp_distribution", "Comp_weights",
                     "comp_weights")$family, 0L)

  # And the diet switch, whose alias is Diet_loglike.
  expect_identical(.rce_switch_column(list(Diet_loglike = 1L),
                                      "Diet_distribution"), 1L)
  expect_null(.rce_switch_column(list(), "Diet_distribution"))
})
