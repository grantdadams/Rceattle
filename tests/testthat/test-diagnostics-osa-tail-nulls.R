# The OSA tail statistic and its null interval must come from the SAME
# estimator. Until v5.29.0 the observed value was stats::quantile()'s type-7
# interpolation while the null was simulated from type-7 quantiles too -- self
# consistent, but seed-dependent. Replacing only the null with a closed form
# would have mixed two estimators: measured coverage of the nominal 95%
# interval falls to about 0.86 near n = 50, and non-monotonically in n, so an
# analyst could not learn to discount it.
#
# Both sides are now the r-th order statistic, whose exact distribution is
# Beta(r, n - r + 1) on the probability scale. The coverage test below is the
# one that would have caught the rejected pairing.

test_that("the order-statistic null interval is the exact Beta result", {
  n <- 200L; probs <- c(0.025, 0.975)
  s <- .osa_tail_null(0.025, n, probs)

  r <- as.integer(round(0.025 * (n + 1)))
  expect_identical(s$r, r)
  expect_equal(s$nominal, r / (n + 1))
  expect_equal(s$null, stats::qnorm(stats::qbeta(probs, r, n - r + 1)))
})


test_that("a short series gets a finite interval rather than a silent pass", {
  # Unclamped, round(0.975 * (n + 1)) exceeds n for n <= 19, and
  # qbeta(p, n + 1, 0) returns 1 -- so qnorm() returns Inf and the upper tail
  # check would pass for every such series. That is a false clean bill of
  # health, which is worse than a missing one.
  for (n in c(3L, 5L, 10L, 19L)) {
    s <- .osa_tail_null(0.975, n, c(0.025, 0.975))
    expect_true(all(is.finite(s$null)), info = paste("n =", n))
    expect_lte(s$r, n)
    expect_gte(s$r, 1L)
  }
  # The clamp binds only where it must: n = 20 is the first unclamped n.
  expect_identical(.osa_tail_null(0.975, 20L)$r, 20L)
  expect_identical(.osa_tail_null(0.975, 200L)$r,
                   as.integer(round(0.975 * 201)))
})


test_that("the tail interval covers at its nominal rate across n", {
  # The rejected type-7 pairing fails at n = 40, 50, 80 and 120; those are kept
  # here as the regression anchors.
  skip_on_cran()
  set.seed(20260909)
  nrep <- 4000

  for (n in c(10L, 25L, 40L, 50L, 80L, 120L, 200L)) {
    lo <- .osa_tail_null(0.025, n)
    hi <- .osa_tail_null(0.975, n)
    cov_lo <- cov_hi <- 0L
    for (i in seq_len(nrep)) {
      s <- sort(stats::rnorm(n))
      if (s[lo$r] >= lo$null[1] && s[lo$r] <= lo$null[2]) cov_lo <- cov_lo + 1L
      if (s[hi$r] >= hi$null[1] && s[hi$r] <= hi$null[2]) cov_hi <- cov_hi + 1L
    }
    expect_gt(cov_lo / nrep, 0.94)
    expect_gt(cov_hi / nrep, 0.94)
    expect_lt(cov_lo / nrep, 0.96)
    expect_lt(cov_hi / nrep, 0.96)
  }
})


test_that("osa_diagnostics() no longer depends on a seed", {
  set.seed(1)
  osa <- data.frame(source = "index", fleet = 1L,
                    residual = stats::rnorm(120))
  a <- suppressWarnings(osa_diagnostics(osa, seed = 1))
  b <- suppressWarnings(osa_diagnostics(osa, seed = 99999))
  expect_identical(as.data.frame(a), as.data.frame(b))
})


test_that("nsim and seed are kept but announced as ignored", {
  osa <- data.frame(source = "index", fleet = 1L, residual = stats::rnorm(60))
  expect_warning(osa_diagnostics(osa, nsim = 500), "ignored")
  expect_warning(osa_diagnostics(osa, seed = 7), "ignored")
  # The default call stays quiet.
  expect_silent(osa_diagnostics(osa))
})


test_that("the reported tail carries the order statistic it actually used", {
  osa <- data.frame(source = "index", fleet = 1L, residual = stats::rnorm(80))
  d <- osa_diagnostics(osa)
  expect_true(all(c("lower_r", "upper_r", "lower_p", "upper_p") %in% names(d)))
  expect_equal(d$lower_p[1], d$lower_r[1] / (d$n[1] + 1))
  # The reported statistic IS that order statistic, not an interpolation.
  expect_equal(d$lower[1], sort(osa$residual)[d$lower_r[1]])
  expect_equal(d$upper[1], sort(osa$residual)[d$upper_r[1]])
})
