# .dsem_draw_projection() -- the conditional GMRF draw that supplies a DSEM's
# projected recruitment deviations in sample_rec().
#
# These test the DISTRIBUTION, not the plumbing, against a precision built by
# hand: given the fitted hindcast states, the projection block must come from
#
#   x_P | x_H  ~  N( mu_P - Q_PP^-1 Q_PH (x_H - mu_H),  Q_PP^-1 )
#
# which is the conditional of the same GMRF that scored the fit. A draw that
# ignored the coupling -- iid, or the marginal rather than the conditional --
# passes no check below. Building Q here rather than fitting keeps this fast and
# makes the expected answer analytic. The end-to-end check -- that the model's
# OWN reported Q and map reach this code, on a real covariate-carrying fit --
# is tools/verify/verify-dsem-projection-draw.R.

# AR1 precision on n states with innovation SD sigma, as an unnormalized GMRF.
.ar1_prec <- function(n, rho, sigma) {
  Q <- matrix(0, n, n)
  for (i in seq_len(n)) Q[i, i] <- if (i %in% c(1L, n)) 1 else 1 + rho^2
  Q[1, 1] <- 1
  for (i in seq_len(n - 1L)) { Q[i, i + 1L] <- -rho; Q[i + 1L, i] <- -rho }
  Q / sigma^2
}

# A minimal object carrying exactly what the draw reads. `fixed` marks nodes the
# map pins (NA in a TMB map), i.e. known values to condition on.
.fake_dsem_fit <- function(n_t, n_j, n_hind, Q, mu, X, fixed = NULL) {
  mp <- matrix(seq_len(n_t * n_j), n_t, n_j)
  if (!is.null(fixed)) mp[fixed] <- NA
  structure(list(
    data_list = list(styr = 1L, endyr = n_hind),
    estimated_params = list(dsem_x_tj = X),
    map = list(mapList = list(dsem_x_tj = mp)),
    quantities = list(dsem_Q = Q, dsem_xhat_tj = mu,
                      dsem_delta_tj = matrix(0, n_t, n_j))),
    class = "Rceattle")
}

testthat::test_that(".dsem_draw_projection() draws the GMRF conditional", {
  n_t <- 30L; n_hind <- 20L; rho <- 0.7; sigma <- 0.5
  Q <- .ar1_prec(n_t, rho, sigma)
  mu <- matrix(0, n_t, 1L)
  X  <- matrix(0, n_t, 1L)
  # A large terminal hindcast state: an AR1 must carry it into the projection
  # and decay it, so this separates the conditional from the marginal.
  X[n_hind, 1L] <- 3
  fit <- .fake_dsem_fit(n_t, 1L, n_hind, Q, mu, X)

  P <- (n_hind + 1L):n_t
  H <- seq_len(n_hind)
  cond_mu <- -solve(Q[P, P], Q[P, H] %*% X[H, 1L])
  cond_sd <- sqrt(diag(solve(Q[P, P])))

  set.seed(42)
  D <- replicate(4000, .dsem_draw_projection(fit)[P, 1L])

  # The mean is the conditional mean, within Monte Carlo error.
  testthat::expect_lt(max(abs(rowMeans(D) - cond_mu) / cond_sd), 4 / sqrt(4000))
  # ... and the spread is the conditional SD, not the marginal one.
  testthat::expect_lt(max(abs(apply(D, 1, stats::sd) / cond_sd - 1)), 0.06)

  # The hindcast is left exactly as fitted -- the draw conditions on it.
  testthat::expect_identical(.dsem_draw_projection(fit)[H, 1L], X[H, 1L])

  # The conditional mean decays the terminal state at rate rho, so the process
  # genuinely propagates rather than reverting to mu immediately.
  testthat::expect_equal(as.numeric(cond_mu[1]), rho * 3, tolerance = 1e-8)
  testthat::expect_gt(cond_mu[1], cond_mu[5])
})

testthat::test_that("sample = FALSE returns the conditional mean, not a draw", {
  n_t <- 24L; n_hind <- 16L
  Q <- .ar1_prec(n_t, 0.6, 0.4)
  mu <- matrix(0.2, n_t, 1L); X <- matrix(0.2, n_t, 1L); X[n_hind, 1L] <- 1.5
  fit <- .fake_dsem_fit(n_t, 1L, n_hind, Q, mu, X)

  a <- .dsem_draw_projection(fit, sample = FALSE)
  b <- .dsem_draw_projection(fit, sample = FALSE)
  testthat::expect_identical(a, b)

  P <- (n_hind + 1L):n_t; H <- seq_len(n_hind)
  cond_mu <- mu[P, 1L] - solve(Q[P, P], Q[P, H] %*% (X[H, 1L] - mu[H, 1L]))
  testthat::expect_equal(a[P, 1L], as.numeric(cond_mu), tolerance = 1e-10)
})

testthat::test_that("a coupled column moves with its partner, an uncoupled one does not", {
  # Two columns, the first AR1 and the second white noise. The conditional must
  # propagate the first column's terminal state and leave the second's mean at
  # mu -- an iid draw would fail the first, and a draw that smoothed everything
  # would fail the second.
  n_t <- 20L; n_hind <- 14L
  Q <- as.matrix(Matrix::bdiag(.ar1_prec(n_t, 0.8, 0.5),
                               diag(1 / 0.5^2, n_t)))
  mu <- matrix(0, n_t, 2L); X <- matrix(0, n_t, 2L)
  X[n_hind, 1L] <- 2; X[n_hind, 2L] <- 2
  fit <- .fake_dsem_fit(n_t, 2L, n_hind, Q, mu, X)

  m <- .dsem_draw_projection(fit, sample = FALSE)
  testthat::expect_equal(m[n_hind + 1L, 1L], 0.8 * 2, tolerance = 1e-8)
  testthat::expect_equal(m[n_hind + 1L, 2L], 0, tolerance = 1e-10)
})

testthat::test_that("a hindcast-only latent field is refused, not extrapolated", {
  n_t <- 12L
  Q <- .ar1_prec(n_t, 0.5, 0.5)
  fit <- .fake_dsem_fit(n_t, 1L, n_t, Q, matrix(0, n_t, 1L), matrix(0, n_t, 1L))
  testthat::expect_error(.dsem_draw_projection(fit),
                         "estimate_projection = TRUE")

  # A precision that does not match the latent field is a wiring error, not
  # something to align by guesswork.
  bad <- .fake_dsem_fit(n_t, 1L, 8L, .ar1_prec(n_t - 1L, 0.5, 0.5),
                        matrix(0, n_t, 1L), matrix(0, n_t, 1L))
  testthat::expect_error(.dsem_draw_projection(bad), "cannot be aligned")
})

testthat::test_that("a known future covariate is conditioned on, not drawn away", {
  # THE case a climate-linked DSEM exists for. Under family = "fixed" a
  # covariate column of x_tj IS the environmental data, pinned by the map. Where
  # env_data runs past endyr those projection nodes are a KNOWN future scenario.
  # Drawing them instead of conditioning on them replaces the scenario with
  # noise and returns the no-covariate answer -- projected recruitment stops
  # responding to the environment at all, silently.
  #
  # Two columns: recdevs (free everywhere) and env (fixed everywhere), with a
  # lag-0 path env -> recdev of coefficient beta. Then
  #   E[recdev_t | env] = beta * env_t
  # exactly, so a correct conditional tracks the covariate and a draw that
  # includes the env nodes in P returns ~0.
  n_t <- 20L; n_hind <- 12L; beta <- 0.75; sd_r <- 0.4; sd_e <- 1
  n <- n_t * 2L
  # Q = (I - Rho)' V^-1 (I - Rho), with the single path env_t -> recdev_t.
  Rho <- matrix(0, n, n)
  for (t in seq_len(n_t)) Rho[t, n_t + t] <- beta
  V <- diag(c(rep(sd_r^2, n_t), rep(sd_e^2, n_t)))
  IR <- diag(n) - Rho
  Q <- t(IR) %*% solve(V) %*% IR

  env <- as.numeric(scale(sin(seq_len(n_t))))
  X <- cbind(rep(0, n_t), env)
  X[seq_len(n_hind), 1L] <- beta * env[seq_len(n_hind)]   # fitted hindcast
  mu <- matrix(0, n_t, 2L)
  # The env column is fixed in EVERY year -- that is what family = "fixed" does.
  fit <- .fake_dsem_fit(n_t, 2L, n_hind, Q, mu, X,
                        fixed = cbind(seq_len(n_t), 2L))

  m <- .dsem_draw_projection(fit, sample = FALSE)
  pj <- (n_hind + 1L):n_t

  # The recruitment states track the known future environment...
  testthat::expect_equal(m[pj, 1L], beta * env[pj], tolerance = 1e-8)
  # ... which is a real signal, not a flat zero the buggy version would give.
  testthat::expect_gt(max(abs(beta * env[pj])), 0.5)
  # ... and the environment itself is left exactly as supplied, never redrawn.
  testthat::expect_identical(m[, 2L], X[, 2L])

  # The same holds for a sampled draw: the env column must not move, in any
  # draw, and the recdevs must centre on the covariate response.
  set.seed(9)
  D <- replicate(400, .dsem_draw_projection(fit)[, 1L][pj])
  E <- replicate(20, .dsem_draw_projection(fit)[, 2L])
  testthat::expect_true(all(apply(E, 2, identical, X[, 2L])))
  testthat::expect_lt(max(abs(rowMeans(D) - beta * env[pj])), 4 * sd_r / sqrt(400))
  # The draw variance is the CONDITIONAL one. Drawing the env column too would
  # inflate it to sd_r^2 + beta^2 * sd_e^2 -- and the template's bias correction
  # is computed conditional on the environment, so that inflation would fight it.
  testthat::expect_lt(max(abs(apply(D, 1, stats::sd) / sd_r - 1)), 0.12)
  testthat::expect_lt(mean(apply(D, 1, stats::var)),
                      sd_r^2 + 0.5 * beta^2 * sd_e^2)
})

testthat::test_that("a fully-fixed projection is refused rather than returned unchanged", {
  n_t <- 16L; n_hind <- 10L
  Q <- .ar1_prec(n_t, 0.5, 0.5)
  fit <- .fake_dsem_fit(n_t, 1L, n_hind, Q, matrix(0, n_t, 1L),
                        matrix(0, n_t, 1L),
                        fixed = cbind((n_hind + 1L):n_t, 1L))
  testthat::expect_error(.dsem_draw_projection(fit), "nothing to draw")
})
