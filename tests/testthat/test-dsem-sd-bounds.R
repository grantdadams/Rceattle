# Bounding a DSEM's standard deviations away from the mirrored optimum.
#
# PROVENANCE. dsem_beta_z holds every free SEM path in one vector: `->` paths
# are regression coefficients, `<->` paths are entries of the Cholesky factor of
# the exogenous covariance. A lag-0 two-headed SELF path is a diagonal of that
# factor -- a standard deviation -- and the likelihood sees only Gamma'Gamma, so
# its sign is not identified; the template reads it as sqrt(square(beta_z)).
# Unbounded, the surface is exactly symmetric about 0. For an MLE that is a
# harmless pair of mirrored optima. For MCMC it is fatal: the posterior is
# bimodal by construction, so chains do not mix, R-hat is meaningless and any
# interval on sigma is wrong.
#
# The restriction is only exact for a diagonal that stands ALONE in its row. A
# cross-covariance (A <-> B) or a lagged two-headed path puts off-diagonal
# entries in that row, sign-flipping one element is no longer invariant, and
# identifying it means flipping the whole row -- which is not attempted, so
# those keep unbounded support.

.sem_full <- function(...) {
  rows <- list(...)
  out <- do.call(rbind, lapply(rows, function(r) as.data.frame(r,
                                                stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}
.row <- function(first, second, direction, lag, parameter) {
  list(path = paste(first, if (direction == 2) "<->" else "->", second),
       lag = lag, name = paste0("p", parameter), start = 0.5,
       parameter = parameter, first = first, second = second,
       direction = direction)
}

testthat::test_that("only standalone lag-0 variance paths are treated as SDs", {
  sf <- .sem_full(
    .row("recdevs1", "recdevs1", 1, 1, 1),   # -> self path: an AR coefficient
    .row("temp",     "recdevs1", 1, 0, 2),   # -> covariate effect
    .row("recdevs1", "recdevs1", 2, 0, 3),   # <-> lag 0 self: an SD
    .row("recdevs2", "recdevs2", 2, 0, 4))   # <-> lag 0 self: an SD
  got <- Rceattle:::.dsem_sd_indices(list(sem_full = sf))
  testthat::expect_identical(got$sd, c(3L, 4L))
  testthat::expect_length(got$entangled, 0L)
})

testthat::test_that("a cross-covariance leaves both variables' SDs unbounded", {
  # Gamma's row for recdevs1 now has an off-diagonal entry, so flipping the sign
  # of its diagonal alone changes Gamma'Gamma -- the surface is not symmetric
  # and the bound would cut off a region the likelihood can distinguish.
  sf <- .sem_full(
    .row("recdevs1", "recdevs1", 2, 0, 1),
    .row("recdevs2", "recdevs2", 2, 0, 2),
    .row("recdevs3", "recdevs3", 2, 0, 3),
    .row("recdevs1", "recdevs2", 2, 0, 4))
  got <- Rceattle:::.dsem_sd_indices(list(sem_full = sf))
  testthat::expect_identical(got$sd, 3L)
  testthat::expect_setequal(got$entangled, c("recdevs1", "recdevs2"))
})

testthat::test_that("a lagged variance path is not a standard deviation", {
  # `x <-> x, 1` is a covariance between successive years -- off-diagonal in the
  # stacked field, whatever the variable names say.
  sf <- .sem_full(
    .row("recdevs1", "recdevs1", 2, 1, 1),
    .row("recdevs2", "recdevs2", 2, 0, 2))
  got <- Rceattle:::.dsem_sd_indices(list(sem_full = sf))
  testthat::expect_identical(got$sd, 2L)
})

testthat::test_that("a parameter shared with a non-variance path is not bounded", {
  # Two paths constrained to one name share a beta_z entry. Bounding it because
  # one of them is a variance would bound the other one too.
  sf <- .sem_full(
    .row("recdevs1", "recdevs1", 2, 0, 1),
    .row("temp",     "recdevs1", 1, 0, 1),
    .row("recdevs2", "recdevs2", 2, 0, 2))
  got <- Rceattle:::.dsem_sd_indices(list(sem_full = sf))
  testthat::expect_identical(got$sd, 2L)
})

testthat::test_that("a fixed path (parameter 0) is not indexed", {
  sf <- .sem_full(
    .row("recdevs1", "recdevs1", 2, 0, 0),
    .row("recdevs2", "recdevs2", 2, 0, 1))
  got <- Rceattle:::.dsem_sd_indices(list(sem_full = sf))
  testthat::expect_identical(got$sd, 1L)
})

testthat::test_that("build_bounds() bounds the SDs and leaves everything else alone", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  idx <- Rceattle:::.dsem_sd_indices(fit$dsem)$sd
  testthat::expect_gt(length(idx), 0L)

  b <- Rceattle::build_bounds(param_list = fit$estimated_params,
                              data_list = fit$data_list, dsem = fit$dsem)
  testthat::expect_true(all(b$lower$dsem_beta_z[idx] == 0))
  testthat::expect_true(all(is.infinite(b$upper$dsem_beta_z)))
  other <- setdiff(seq_along(b$lower$dsem_beta_z), idx)
  testthat::expect_true(all(is.infinite(b$lower$dsem_beta_z[other])))

  # `dsem` is optional, so an existing caller gets exactly what it always did.
  b0 <- Rceattle::build_bounds(param_list = fit$estimated_params,
                               data_list = fit$data_list)
  testthat::expect_true(all(is.infinite(b0$lower$dsem_beta_z)))
})

testthat::test_that("a fit warm-started at the negative root is flipped, not rejected", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # A model fitted before the bound existed is as likely to sit at -sigma as
  # +sigma. Erroring on "inits not within bounds" would strand those fits; the
  # flip is exact where the bound applies, so the refit finds the same optimum.
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  testthat::expect_true(all(as.numeric(fit$estimated_params$dsem_beta_z) >= 0))

  neg <- fit$estimated_params
  neg$dsem_beta_z <- -abs(neg$dsem_beta_z)
  again <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = neg, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  testthat::expect_equal(again$opt$objective, fit$opt$objective, tolerance = 1e-4)
  testthat::expect_true(all(as.numeric(again$estimated_params$dsem_beta_z) >= 0))
})

testthat::test_that("a build_DSEM() spec passed to build_bounds() says it is the wrong object", {
  # Only the BUILT dsem a fit carries has sem_full. Handing in the
  # specification instead is easy to do and, left silent, leaves an MCMC
  # bimodal with nothing pointing at the cause.
  pl <- list(dsem_beta_z = c(0.3, 0.5), log_gam_a = 1, log_gam_b = 1)
  dl <- list(nspp = 1L, suitMode = 0L)
  testthat::expect_warning(
    b <- Rceattle::build_bounds(param_list = pl, data_list = dl,
                                dsem = Rceattle::build_DSEM()),
    "BUILT DSEM")
  testthat::expect_true(all(is.infinite(b$lower$dsem_beta_z)))
})
