# The 2D-AR1 densities, pinned against an independent multivariate-normal
# reference rather than against themselves.
#
# TMB's SEPARABLE(f, g) applies f to the OUTERMOST array dimension and g to the
# fastest-running one (density.hpp). Both fields here are laid out (age, year)
# or (bin, year), so the YEAR correlation has to be passed first. Written the
# other way round the two correlations are silently exchanged: the fit stays
# self-consistent, so nothing errors and no golden model notices -- the reported
# age correlation is simply the year correlation. That is what these pin.
#
# `calc_nll_ar1_2d()` (helpers.R) builds the covariance as
# kronecker(R_year, R_age), i.e. rows = age, so it is the independent statement
# of the intended orientation.

testthat::test_that("the 2D-AR1 M density correlates ages by rho_age and years by rho_year", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mvtnorm")

  d <- make_test_data()
  fit <- suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0),
    M1Fun = Rceattle::build_M1(M1_model = 1, M1_re = 6)))

  SIG <- 0.3; RHO_A <- 0.197; RHO_Y <- 0.862   # deliberately unequal
  nages  <- d$nages[1]
  nyrs   <- d$endyr - d$styr + 1L

  obj <- fit$obj
  par <- obj$env$last.par.best
  if (is.null(par)) { obj$fn(obj$par); par <- obj$env$last.par.best }
  par[names(par) == "M1_dev_log_sd"] <- log(SIG)
  # M1_rho slot 1 = age, slot 2 = year (the 1-D M modes pin this reading:
  # the by-age mode uses slot 1, the by-year mode slot 2).
  par[names(par) == "M1_rho"] <- c(atanh(RHO_A), atanh(RHO_Y))

  # A known deviation field, so the comparison is one number against one number.
  # Written straight into the parameter vector and read back with obj$report(),
  # which evaluates at exactly these values -- obj$fn() would run the Laplace
  # inner problem and move the random effects off the field under test.
  set.seed(11)
  field <- matrix(stats::rnorm(nages * nyrs), nages, nyrs)
  dev_pos <- which(names(par) == "log_M1_dev")
  arr <- array(0, dim = dim(fit$estimated_params$log_M1_dev))
  arr[1, 1, seq_len(nages), seq_len(nyrs)] <- field
  # log_M1_dev is mapped, so only the free cells appear in par; take them in the
  # same order the map lays them out.
  mp <- as.integer(fit$map$mapList$log_M1_dev)
  par[dev_pos] <- vapply(seq_along(dev_pos), function(k)
    arr[which(mp == k)[1]], numeric(1))

  # jnll_comp is the raw reported matrix here: no rownames (rename_output()
  # adds those later), so address the row by its JnllRow constant.
  JNLL_M_RE <- 15L
  got <- obj$report(par)$jnll_comp[JNLL_M_RE + 1L, 1]

  want  <- calc_nll_ar1_2d(field, SIG, rho_a = RHO_A, rho_y = RHO_Y)
  wrong <- calc_nll_ar1_2d(field, SIG, rho_a = RHO_Y, rho_y = RHO_A)

  testthat::expect_equal(got, want, tolerance = 1e-8)
  # The check only means something if the two orientations actually differ.
  testthat::expect_gt(abs(want - wrong), 1)
})
