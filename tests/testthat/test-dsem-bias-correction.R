# The lognormal bias correction under a DSEM: rec_dev = x_tj - bias * Var/2,
# where Var is the variance of the recruitment-deviation state GIVEN whatever
# environmental values the model was handed (ceattle.cpp section 5.5b,
# calculate_dsem()'s margvar block).
#
# "Given whatever it was handed" is per CELL. build_dsem_objects() pads env_data
# to styr:dsem_endyr, so a covariate that starts late, has a gap, or carries no
# projection scenario is a free latent state in those years, and a correction
# that treated the whole column as known would under-correct exactly there --
# measured at 17.1% high on projected recruitment before this was fixed. That is
# invisible to a covariate-free sem, which is why these carry a covariate.
#
# The expected value is read off the model's OWN reported precision:
# Var(x_U | x_K) = diag(solve(Q[U, U])) with U the cells the map leaves free.
# Q and margvar_tj come from different code paths inside calculate_dsem(), so
# agreement is a real cross-check rather than a restatement.

.bias_fit <- function(env, sem) {
  d <- Rceattle::BS2017SS
  d$projyr <- d$endyr + 8L
  d$env_data <- data.frame(Year = d$styr:d$projyr, BT = env)
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0,
    recFun = Rceattle::build_srr(proj_mean_rec = FALSE),
    dsem = Rceattle::build_DSEM(sem = sem, family = "fixed",
                                estimate_projection = TRUE),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
}

# Var(x_U | x_K) at this species' recruitment cells, from the reported precision.
.cond_var_rec <- function(fit, sp = 1L) {
  Q    <- as.matrix(fit$quantities$dsem_Q)
  X    <- fit$estimated_params$dsem_x_tj
  n_t  <- nrow(X); n_j <- ncol(X)
  rcol <- fit$dsem$tmb_inputs$data$rec_dev_col[sp] + 1L
  free <- which(!is.na(as.numeric(fit$map$mapList$dsem_x_tj)))
  S    <- solve(Q[free, free, drop = FALSE])
  vapply(seq_len(n_t), function(t) {
    i <- match((rcol - 1L) * n_t + t, free)
    if (is.na(i)) NA_real_ else S[i, i]
  }, numeric(1))
}

.rec_margvar <- function(fit, sp = 1L) {
  rcol <- fit$dsem$tmb_inputs$data$rec_dev_col[sp] + 1L
  as.matrix(fit$quantities$dsem_margvar_tj)[, rcol]
}

# rho > 0 and a covariate effect, or the conditional and the joint variance
# coincide and nothing below can fail.
SEM_PLAIN <- "recdevs1 <-> recdevs1, 0, sigmaR1, 0.55
recdevs2 <-> recdevs2, 0, sigmaR2, 0.55
recdevs3 <-> recdevs3, 0, sigmaR3, 0.55
recdevs1 -> recdevs1, 1, rho1, 0.6
BT <-> BT, 0, sdBT, 1
BT -> recdevs1, 0, bBT, 0.45"

# The same, with the covariate autoregressive. This is what forces the Schur
# complement: an unobserved BT year drives an observed one, so the conditional
# variance is not "drop the known cells' innovations from the joint".
SEM_AR <- paste0(SEM_PLAIN, "\nBT -> BT, 1, arBT, 0.5")


testthat::test_that("the correction uses the variance given the observed covariate", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d0 <- Rceattle::BS2017SS
  yrs <- d0$styr:(d0$endyr + 8L)
  env <- as.numeric(scale(sin(seq_along(yrs) / 3)))

  fit <- .bias_fit(env, SEM_PLAIN)
  testthat::expect_equal(.rec_margvar(fit), .cond_var_rec(fit), tolerance = 1e-10)

  # With BT known in every year the correction is the conditional variance,
  # sigmaR^2/(1-rho^2), and carries no covariate term at all.
  b   <- fit$estimated_params$dsem_beta_z
  nm  <- fit$dsem$sem_full$name
  idx <- function(x) fit$dsem$sem_full$parameter[nm == x]
  mid <- round(length(yrs) / 2)
  testthat::expect_equal(.rec_margvar(fit)[mid],
                         b[idx("sigmaR1")]^2 / (1 - b[idx("rho1")]^2),
                         tolerance = 1e-8)
})


testthat::test_that("a covariate is known only in the years it was supplied", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d0  <- Rceattle::BS2017SS
  yrs <- d0$styr:(d0$endyr + 8L)
  env <- as.numeric(scale(sin(seq_along(yrs) / 3)))
  # Unobserved at both ends: a late-starting series with no projection scenario,
  # which is the ordinary shape of an ESP covariate.
  gap <- env
  gap[yrs < d0$styr + 15L | yrs > d0$endyr] <- NA_real_

  fit <- .bias_fit(gap, SEM_AR)
  mv  <- .rec_margvar(fit)
  testthat::expect_equal(mv, .cond_var_rec(fit), tolerance = 1e-10)

  # ... and the test has teeth. Where BT is known the variance settles at
  # sigmaR^2/(1-rho^2), which is what a per-column rule returned EVERYWHERE.
  # Where BT is not known it must be strictly larger, because the covariate is
  # then part of the spread rather than part of the mean. Checked away from the
  # known/unknown boundaries, where the recruitment AR1 carries variance across
  # and neither value is stationary.
  n_t   <- length(yrs)
  b   <- fit$estimated_params$dsem_beta_z
  nm  <- fit$dsem$sem_full$name
  idx <- function(x) fit$dsem$sem_full$parameter[nm == x]
  per_column <- b[idx("sigmaR1")]^2 / (1 - b[idx("rho1")]^2)

  mid_known <- which(yrs == d0$endyr - 5L)
  testthat::expect_equal(mv[mid_known], per_column, tolerance = 1e-6)
  # The terminal projection year -- the one that sets an ABC -- and a hindcast
  # year before the covariate starts.
  testthat::expect_gt(mv[n_t] / per_column, 1.05)
  testthat::expect_gt(mv[which(yrs == d0$styr + 10L)] / per_column, 1.05)
})


testthat::test_that("a cell the model was given carries no correction", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d0  <- Rceattle::BS2017SS
  yrs <- d0$styr:(d0$endyr + 8L)
  env <- as.numeric(scale(sin(seq_along(yrs) / 3)))

  fit  <- .bias_fit(env, SEM_PLAIN)
  mv   <- as.matrix(fit$quantities$dsem_margvar_tj)
  n_t  <- nrow(mv)
  rcol <- fit$dsem$tmb_inputs$data$rec_dev_col[1] + 1L
  pinned <- colSums(matrix(is.na(as.numeric(fit$map$mapList$dsem_x_tj)),
                           n_t, ncol(mv)))
  ecol <- setdiff(which(pinned > 0), rcol)

  # BT is data in every year, so it has no variance left to report.
  testthat::expect_length(ecol, 1L)
  testthat::expect_true(all(mv[, ecol] == 0))
  # The recruitment columns are latent throughout, so they must not be zeroed.
  testthat::expect_true(all(mv[, rcol] > 0))
})
