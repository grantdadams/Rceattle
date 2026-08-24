# process_residuals(process = "recruitment") on a DSEM.
#
# A DSEM's recruitment deviations are the latent states of a GMRF, so there is
# no per-year normal prior to divide by. The residual is defined all the same:
# whiten the states with the process precision. Conditioning on the states the
# model was GIVEN, the estimated states have precision Q[E, E] and mean
# mu_E - Q[E,E]^-1 Q[E,F] (x_F - mu_F), and r = chol(Q[E,E]) %*% (x_E - m) has
# covariance I -- so under a correctly specified sem the residuals are
# approximately iid N(0, 1), the same claim the scalar path makes.
#
# Two things are pinned here. The ALGEBRA, against draws from the process prior,
# which needs no data and no fit and would catch a transposed Cholesky, a
# mis-signed conditioning term or a column-major indexing slip. And the CELL
# SELECTION: dsem does NOT map an observed covariate cell off, so selecting by
# the parameter map returns one residual per unobserved covariate year as well,
# labelled "recruitment" and pooled by osa_diagnostics(). Measured on a 49-year
# GOApollock fit with one partially observed covariate: 58 rows instead of 49.

testthat::test_that("the whitening algebra returns iid standard normals", {
  # A small GMRF standing in for the DSEM's precision: an AR1 over 30 steps
  # times two variables, so the conditioning block is non-empty and the
  # precision is genuinely correlated.
  n_t <- 30L; n_j <- 2L; n_k <- n_t * n_j
  rho <- 0.6; sig <- 0.8
  A <- diag(n_t); for (t in 2:n_t) A[t, t - 1] <- -rho
  Q1 <- crossprod(A) / sig^2
  Q  <- as.matrix(Matrix::bdiag(Q1, Q1))
  mu <- rep(0, n_k)

  E <- 1:n_t                 # "recruitment" column
  F_ <- setdiff(seq_len(n_k), E)
  xF <- stats::rnorm(length(F_))
  m  <- mu[E] - solve(Q[E, E], Q[E, F_] %*% (xF - mu[F_]))
  U  <- chol(Q[E, E])

  set.seed(11)
  R <- replicate(3000, {
    x <- as.numeric(m) + backsolve(U, stats::rnorm(length(E)))
    as.numeric(U %*% (x - m))
  })
  testthat::expect_equal(mean(R), 0, tolerance = 0.03)
  testthat::expect_equal(stats::sd(as.numeric(R)), 1, tolerance = 0.03)
  # Every element standard, and no correlation left between them.
  testthat::expect_true(all(abs(apply(R, 1, stats::sd) - 1) < 0.1))
  cc <- stats::cor(t(R)); diag(cc) <- 0
  testthat::expect_lt(max(abs(cc)), 0.12)
})

testthat::test_that("a DSEM's recruitment residuals are one per hindcast year", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not_installed("TMB")

  d <- Rceattle::GOApollock
  d$projyr <- 2020
  # Partially observed on purpose: those covariate years are latent states the
  # map does not pin, and are exactly what a map-based selection would return
  # as extra "recruitment" residuals.
  d$env_data$ScaledBT <- as.numeric(scale(d$env_data$BTempC))
  sem <- "ScaledBT -> ScaledBT, 1, AR_BT, 0
ScaledBT -> recdevs1, 1, BT_to_R, 0
recdevs1 <-> recdevs1, 0, sigmaR1, 1"
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, initMode = 2,
    dsem = Rceattle::build_DSEM(sem = sem, family = "fixed"),
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0,
                                        getsd = TRUE))))

  r <- suppressWarnings(Rceattle::process_residuals(fit, process = "recruitment"))
  nyrs <- d$endyr - d$styr + 1L
  testthat::expect_equal(nrow(r), nyrs)
  testthat::expect_setequal(r$year, d$styr:d$endyr)
  testthat::expect_true(all(r$species == 1L))
  testthat::expect_true(all(is.finite(r$residual)))
  testthat::expect_true(inherits(r, "rceattle_osa"))

  # "all" reaches every process now, rather than refusing because one of them
  # cannot be computed.
  a <- suppressWarnings(Rceattle::process_residuals(fit, process = "all"))
  testthat::expect_setequal(unique(a$source),
                            c("recruitment", "initial", "catchability"))

  # "initial" standardizes by the SD the density used, which under a DSEM comes
  # from the sem's variance path -- not the mapped-out exp(R_log_sd) placeholder.
  i <- suppressWarnings(Rceattle::process_residuals(fit, process = "initial"))
  R_sd <- as.numeric(fit$quantities$R_sd)[1]
  R_log_sd <- exp(as.numeric(fit$obj$env$parList()$R_log_sd))[1]
  testthat::expect_false(isTRUE(all.equal(R_sd, R_log_sd)))
  testthat::expect_true(all(is.finite(i$residual)))
})
