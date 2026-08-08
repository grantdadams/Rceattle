# The DSEM recruitment hook, pinned by an exact equivalence.
#
# With a base sem -- only the two-headed self-loop on each species' recdevs --
# the GMRF reduces to N(0, sigma) on the recruitment deviations. That is exactly
# the standard recruitment-deviation density once the lognormal bias correction
# is off (the standard one is N(-sigma^2/2, sigma)). So a base DSEM and a
# non-DSEM fit with bias_adjust_proc = FALSE must agree to the last digit.
#
# This is the strongest available check on the hook: it exercises the
# calculate_dsem() call, the rec_dev overwrite from the latent states, the R_sd
# derivation from beta_z, the JNLL_REC_DEV gate and the JNLL_DSEM accumulation,
# and it fails if any one of them is off by anything at all.

testthat::test_that("a base DSEM equals a non-DSEM fit with bias correction off", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  # Bias correction OFF, so the two share a parameterization: the DSEM's latent
  # state x_tj IS the recruitment deviation. With it ON the model is still the
  # same (rec_dev = x_tj - sigma^2/2, so rec_dev keeps its N(-sigma^2/2, sigma)
  # marginal) but the two are offset in parameter space, so equal parameter
  # values are not equal deviations and a fixed-point comparison is meaningless.
  # The offset itself is checked in the next test.
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0,
                              bias_adjust_proc = FALSE)

  # sem = NULL synthesizes the per-species IID self-loop -- the base DSEM.
  with_dsem <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  without <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, fit_control = fc)))

  a <- with_dsem$quantities
  b <- without$quantities

  # The whole objective, and the states it implies, must be identical.
  testthat::expect_equal(sum(a$jnll_comp), sum(b$jnll_comp), tolerance = 1e-10)
  # Named for what the report actually carries -- `rec_dev` and `biomassSSB` are
  # NOT reported, so asserting on them compares NULL to NULL and passes for free.
  testthat::expect_false(is.null(a$ssb))
  testthat::expect_equal(a$ssb,   b$ssb,   tolerance = 1e-12)
  testthat::expect_equal(a$R,     b$R,     tolerance = 1e-12)
  testthat::expect_equal(a$R_sd,  b$R_sd,  tolerance = 1e-12)

  # ...and the density must have MOVED, not been duplicated or dropped: out of
  # the recruitment-deviation row and into the DSEM row, same value.
  # JNLL_REC_DEV = 10 and JNLL_DSEM = 21 in ceattle.cpp, so rows 11 and 22 here.
  testthat::expect_equal(sum(a$jnll_comp[11L, ]), 0)
  testthat::expect_equal(sum(b$jnll_comp[22L, ]), 0)
  testthat::expect_equal(sum(a$jnll_comp[22L, ]), sum(b$jnll_comp[11L, ]),
                         tolerance = 1e-10)
  testthat::expect_gt(sum(a$jnll_comp[22L, ]), 0)   # and it is not trivially zero
})

testthat::test_that("bias correction is applied to the deviations under a DSEM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # Without this the GMRF would centre the deviations at 0 while the standard
  # density centres them at -sigma^2/2 -- so giving a model a DSEM would change
  # its shrinkage target, and its recruitment, SSB and ABC, silently. init_dev
  # keeps its correction either way, so the two must not disagree.
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)

  on  <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3, random_rec = TRUE,
    msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  off <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3, random_rec = TRUE,
    msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0,
                                        bias_adjust_proc = FALSE))))

  # At x_tj = 0 the deviations differ by exactly -sigma^2/2 per species.
  sd_on <- on$quantities$R_sd
  testthat::expect_equal(as.numeric(on$quantities$R), 
                         as.numeric(off$quantities$R) * exp(-sd_on^2 / 2)[1],
                         tolerance = 1e-6)
})

testthat::test_that("the density is exercised away from x_tj = 0, and never double counted", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # estimateMode 3 builds the DSEM in debug mode, which maps every DSEM
  # parameter out and leaves x_tj identically 0 -- there the equivalence above
  # tests only the normalizing constant. Fit properly so the latent states move,
  # then confirm the density is real and still lands in exactly one row.
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  x <- fit$estimated_params$dsem_x_tj
  testthat::expect_gt(length(x), 0)                 # estimated, not mapped out
  testthat::expect_gt(max(abs(x)), 1e-3)            # and it moved off the origin

  testthat::expect_gt(sum(fit$quantities$jnll_comp[22L, ]), 0)
  testthat::expect_equal(sum(fit$quantities$jnll_comp[11L, ]), 0)
  testthat::expect_true(is.finite(fit$opt$objective))
})

testthat::test_that("the DSEM row is exactly zero for every model without one", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))

  testthat::expect_equal(nrow(fit$quantities$jnll_comp), 22L)
  testthat::expect_identical(sum(fit$quantities$jnll_comp[22L, ]), 0)
  testthat::expect_equal(rownames(fit$quantities$jnll_comp)[22L], "DSEM")
})

testthat::test_that("rec_dev and R_log_sd are mapped out under a DSEM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # The template overwrites both from the latent states, so left estimated they
  # would be free parameters with exactly zero gradient -- flat directions and a
  # singular Hessian rather than an error.
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))

  testthat::expect_true(all(is.na(fit$map$mapList$rec_dev)))
  testthat::expect_true(all(is.na(fit$map$mapList$R_log_sd)))
  testthat::expect_false(any(names(fit$obj$par) %in% c("rec_dev", "R_log_sd")))
})
