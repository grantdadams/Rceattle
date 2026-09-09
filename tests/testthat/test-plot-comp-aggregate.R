# The aggregated composition pools counts across years rather than averaging
# proportions, so the aggregate is a sum of independent draws and its variance
# is the exact sum of the per-year variances. afscOSA's qbinom band at the
# pooled sample size was deliberately not ported: it assumes a distribution the
# old mean-of-proportions statistic did not have. Checked by simulation.

set.seed(90210)

# One panel's worth of long-format rows, as .comp_resid_long() produces.
make_long <- function(obs, p, N, sd = NULL) {
  Y <- nrow(p); B <- ncol(p)
  data.frame(
    Fleet = 1L, Fleet_name = "f", Species = 1L, Sex = 0L, comp_type = 0L,
    Year = rep(seq_len(Y), times = B), N = rep(N, times = B),
    bin = rep(seq_len(B), each = Y),
    obs = as.numeric(obs), hat = as.numeric(p),
    Sd = if (is.null(sd)) as.numeric(sqrt(p * (1 - p) / N)) else as.numeric(sd),
    sex_grp = "combined", bin_lab = "Age", panel = "f", type_lab = "age",
    stringsAsFactors = FALSE)
}


test_that("the aggregate pools counts rather than averaging proportions", {
  # Two years with very different sample sizes and different compositions. The
  # pooled statistic must follow the big year; a mean of proportions would sit
  # halfway between them.
  p   <- rbind(c(0.8, 0.1, 0.1), c(0.1, 0.1, 0.8))
  obs <- p
  N   <- c(10, 990)

  agg <- .comp_aggregate(make_long(obs, p, N))
  expect_equal(agg$obs[agg$bin == 1], (0.8 * 10 + 0.1 * 990) / 1000,
               tolerance = 1e-10)
  expect_equal(agg$ISS[1], 1000)
  # A mean of proportions would have given 0.45 in bin 1.
  expect_false(isTRUE(all.equal(agg$obs[agg$bin == 1], 0.45, tolerance = 1e-6)))
})


test_that("the aggregate band covers at its nominal rate", {
  skip_on_cran()
  Y <- 25L; B <- 8L; nrep <- 1500
  p <- t(replicate(Y, { v <- exp(-0.5 * ((1:B - 4) / 1.5)^2) + 0.02; v / sum(v) }))
  N <- round(exp(stats::rnorm(Y, log(150), 0.6)))

  inside <- matrix(NA, nrep, B)
  for (i in seq_len(nrep)) {
    obs <- t(vapply(seq_len(Y),
                    function(y) stats::rmultinom(1, N[y], p[y, ])[, 1] / N[y],
                    numeric(B)))
    a <- .comp_aggregate(make_long(obs, p, N))
    a <- a[order(a$bin), ]
    inside[i, ] <- a$obs >= a$lwr & a$obs <= a$upr
  }
  cover <- colMeans(inside)
  # Normal approximation on an exact variance: close to nominal, and never
  # wildly conservative the way a pooled-binomial band would be.
  expect_gt(min(cover), 0.90)
  expect_lt(max(cover), 0.99)
})


test_that("projection rows are excluded from the aggregate and its ISS", {
  p   <- rbind(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5))
  obs <- p
  N   <- c(100, 100, 100)
  long <- make_long(obs, p, N)
  long$Year <- rep(c(1, 2, 99), times = 2)   # year 99 is past endyr

  expect_equal(.comp_aggregate(long, endyr = 2)$ISS[1], 200)
  expect_equal(.comp_aggregate(long)$ISS[1], 300)   # ungated, for comparison
})


test_that("a bin with no assumed sd keeps its composition, losing only the band", {
  # aggregate()'s formula interface drops whole rows on any NA, so this bin
  # would otherwise vanish from the aggregated composition entirely.
  p   <- rbind(c(0.5, 0.3, 0.2), c(0.5, 0.3, 0.2))
  obs <- p
  N   <- c(100, 100)
  long <- make_long(obs, p, N)
  long$Sd[long$bin == 2] <- NA_real_

  agg <- .comp_aggregate(long)
  expect_equal(sort(agg$bin), 1:3)
  expect_equal(agg$obs[agg$bin == 2], 0.3, tolerance = 1e-10)
  expect_true(is.na(agg$lwr[agg$bin == 2]))
  expect_true(all(is.finite(agg$lwr[agg$bin != 2])))
})


test_that("a wider assumed sd widens the band proportionally", {
  # The band reads its variance from the same Sd the residuals were divided by,
  # so an overdispersed fleet gets a wider band without any change here.
  p   <- rbind(c(0.4, 0.6), c(0.4, 0.6))
  obs <- p
  N   <- c(100, 100)
  narrow <- .comp_aggregate(make_long(obs, p, N))
  wide   <- .comp_aggregate(make_long(obs, p, N,
                                      sd = sqrt(p * (1 - p) / N) * 3))
  expect_equal((wide$upr - wide$lwr) / (narrow$upr - narrow$lwr), rep(3, 2),
               tolerance = 1e-8)
})
