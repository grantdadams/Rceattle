# A NAIVE (IID) DSEM must reproduce a non-DSEM model -- not just in the fit, but
# through the diagnostics.
#
# A base DSEM is only the per-species two-headed self-loop, so the GMRF reduces
# to N(0, sigma) on the recruitment deviations: the standard density once the
# lognormal bias correction is off. The FIT equivalence was already pinned
# (test-dsem-recruitment-hook.R). This pins the RETROSPECTIVE, which was not, and
# which is where two real defects hid:
#
#   1. fit_mod() dropped the DSEM blocks out of a supplied `inits`
#      (build_params() does not declare them, so start_par[names(.skel)] removed
#      them and merge_dsem_params() refilled from the template). Invisible
#      wherever the DSEM is estimated -- any start reaches the same MLE -- and
#      wrong wherever it is PINNED. retrospective()'s forecast refit pins the
#      whole DSEM, so every peel ran with beta_z at its start value: sigma_R =
#      0.707107 for every species, an 80% error in the peeled hindcast, and
#      Mohn's rho flipped sign (-0.49 against +0.13).
#   2. The forecast bias adjustment wrote rec_dev, which the template overwrites
#      from the latent states -- a silent no-op that left the GMRF's own
#      projection in place and moved forecast recruitment ~44%.
#
# Both only show up against the non-DSEM path, which is why comparing a DSEM peel
# to a DSEM refit (as the earlier verification did) could not catch either.

testthat::test_that("a naive DSEM reproduces the non-DSEM retrospective", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)

  # bias_adjust_proc = FALSE on BOTH: the GMRF centres the deviations at 0 while
  # the standard density centres them at -sigma^2/2, so with it on the two are
  # not the same model and the comparison is meaningless.
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0,
                              bias_adjust_proc = FALSE)
  fit_d <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  fit_p <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, fit_control = fc)))

  # Precondition: the fits themselves agree, or a peel difference says nothing.
  testthat::expect_equal(fit_d$opt$objective, fit_p$opt$objective,
                         tolerance = 1e-6)

  rd <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit_d, peels = 1, cores = 1, getsd = FALSE)))
  rp <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit_p, peels = 1, cores = 1, getsd = FALSE)))
  testthat::expect_equal(length(rd$Rceattle_list), length(rp$Rceattle_list))

  pd <- rd$Rceattle_list[[1]]
  pp <- rp$Rceattle_list[[1]]

  # The recruitment SD must be ESTIMATED in the peel, not sitting at a start
  # value. This is the assertion that fails on defect (1): every species came
  # back at exactly 0.707107, identical to six decimals, which no estimator does
  # -- so the spread across species is the sharper check than the values.
  testthat::expect_gt(stats::sd(as.numeric(pd$quantities$R_sd)), 1e-6)

  # The two peels agree on sigma_R to a few percent, not exactly, and the gap is
  # real rather than numerical. The recruitment-deviation density runs over every
  # hindcast year -- endyr_peel appears nowhere in ceattle.cpp -- so a non-DSEM
  # peel pins the peeled-year rec_dev at 0 and STILL SCORES them, contributing
  # ~log(sigma) per peeled year and pulling sigma_R down. A DSEM peel
  # marginalizes those states instead, so it carries no such term. Measured at
  # one peel of 39 years: DSEM 1.0070/1.3465/0.8113 against plain
  # 0.9955/1.3260/0.7990 -- the DSEM value larger in every species, ~1.5%, in the
  # direction the argument predicts. It grows with peel depth, and it is the
  # non-DSEM path that is biased. Everything downstream (hindcast, forecast,
  # Mohn's rho) still agrees to the tolerances below.
  testthat::expect_equal(as.numeric(pd$quantities$R_sd),
                         as.numeric(pp$quantities$R_sd), tolerance = 0.05)
  testthat::expect_true(all(as.numeric(pd$quantities$R_sd) >=
                              as.numeric(pp$quantities$R_sd) - 1e-8))

  keep <- seq_len(pd$data_list$endyr_peel - d$styr + 1L)
  fore <- (length(keep) + 1L):ncol(pd$quantities$R)
  for (sp in seq_len(d$nspp)) {
    # Hindcast: defect (1).
    testthat::expect_equal(as.numeric(pd$quantities$R[sp, keep]),
                           as.numeric(pp$quantities$R[sp, keep]),
                           tolerance = 1e-2, info = paste("hindcast R sp", sp))
    testthat::expect_equal(as.numeric(pd$quantities$ssb[sp, keep]),
                           as.numeric(pp$quantities$ssb[sp, keep]),
                           tolerance = 1e-2, info = paste("hindcast SSB sp", sp))
    # Forecast: defect (2).
    testthat::expect_equal(as.numeric(pd$quantities$R[sp, fore]),
                           as.numeric(pp$quantities$R[sp, fore]),
                           tolerance = 1e-2, info = paste("forecast R sp", sp))
  }

  # NOTE: this is a PEEL-1 contract. The sigma_R gap above grows with peel depth
  # as sqrt((N-k)/N) -- at N = 39 that is -1.3% at k=1, -6.6% at k=5 (the
  # retrospective() default) and -13.8% at k=10 -- so raising `peels` here
  # without widening the tolerances will look like a regression when it is the
  # known non-DSEM bias compounding.

  # And the number a user actually reads.
  num <- vapply(rd$mohns, is.numeric, logical(1))
  testthat::expect_lt(
    max(abs(as.matrix(rd$mohns[, num]) - as.matrix(rp$mohns[, num])), na.rm = TRUE),
    0.05)
})
