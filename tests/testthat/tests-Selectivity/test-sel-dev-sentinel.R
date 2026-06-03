testthat::test_that("Time_varying_sel_sd_prior <= 0 sentinel skips the sel-dev prior", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Setup: IID time-varying DoubleLogistic on a single fleet with non-zero
  # injected sel_inf_dev / log_sel_slp_dev values. With a positive
  # Time_varying_sel_sd_prior, the cpp adds an N(0, sigma) prior on each
  # deviate -> Selectivity_deviates jnll > 0. With the sentinel (<= 0 or NA),
  # build_params sets sel_dev_log_sd = -999 and the cpp gate
  # (sel_dev_log_sd > -100) skips the penalty -> jnll = 0 EXACTLY.
  #
  # The sentinel gates the cpp IID branch (flt_varying_sel == 1 || 2) AND
  # the RandomWalk branch (4 || 5) for sel types in the logistic family
  # (Logistic=1, DoubleLogistic=3, DescendingLogistic=4, DoubleNormal=8).
  # It does NOT apply to NonParametric (2), Hake (5), or 2D/3D AR1 (6/7),
  # which use different penalty structures (Sel_curve_pen / AR1 process).

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 &
                    GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <- 0
  GOA2018SS$fleet_control$Selectivity   <- "DoubleLogistic"
  GOA2018SS$fleet_control$Selectivity_index <- seq_len(nrow(GOA2018SS$fleet_control))
  GOA2018SS$fleet_control$Time_varying_sel  <- "IID"

  # Build params skeleton + inject non-zero sel devs so the prior would fire
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1.0
  mod0  <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params
  # Inject 0.5 into every dev slot (deliberate non-zero so prior would fire)
  inits$sel_inf_dev[]     <- 0.5
  inits$log_sel_slp_dev[] <- 0.5

  # ---- Run A: prior ON (sd = 1.0) ----
  modA <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0)))
  # Selectivity deviates slot (jnll_comp row index 6 in cpp = R row 6)
  rl <- rownames(modA$quantities$jnll_comp)
  sel_dev_row <- grep("Selectivity deviates", rl)
  nll_on <- sum(modA$quantities$jnll_comp[sel_dev_row, ])
  testthat::expect_gt(nll_on, 0)

  # ---- Run B: sentinel (sd = -1) skips the prior ----
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- -1
  mod0B  <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0)))
  initsB <- mod0B$estimated_params
  initsB$sel_inf_dev[]     <- 0.5
  initsB$log_sel_slp_dev[] <- 0.5

  modB <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = initsB, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0)))
  nll_off <- sum(modB$quantities$jnll_comp[sel_dev_row, ])
  testthat::expect_equal(nll_off, 0)
})


testthat::test_that("NA Time_varying_sel_sd_prior also triggers the sentinel", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 &
                    GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <- 0
  GOA2018SS$fleet_control$Selectivity   <- "DoubleLogistic"
  GOA2018SS$fleet_control$Selectivity_index <- seq_len(nrow(GOA2018SS$fleet_control))
  GOA2018SS$fleet_control$Time_varying_sel  <- "IID"
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- NA_real_

  mod0  <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params
  inits$sel_inf_dev[]     <- 0.5
  inits$log_sel_slp_dev[] <- 0.5

  mod <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0)))
  rl <- rownames(mod$quantities$jnll_comp)
  sel_dev_row <- grep("Selectivity deviates", rl)
  testthat::expect_equal(sum(mod$quantities$jnll_comp[sel_dev_row, ]), 0)
})
