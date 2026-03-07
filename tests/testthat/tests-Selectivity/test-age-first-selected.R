testthat::test_that("Test age-based non-parametric selectivity bin-first selected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Adjust data
  n_sel_bins <- 8
  GOA2018SS$fleet_control$Selectivity <- 2
  GOA2018SS$fleet_control$N_sel_bins <- n_sel_bins
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_curve_pen1 <- 5
  GOA2018SS$fleet_control$Sel_curve_pen2 <- 10
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA
  GOA2018SS$fleet_control$Bin_first_selected <- 3
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 3


  # Set params
  inits <- suppressMessages(build_params(GOA2018SS))
  log_selcoffs <- rnorm(n_sel_bins)
  log_selcoffs2 <- rnorm(n_sel_bins)
  inits$sel_coff[,1,1:8] <- rep(log_selcoffs, each = dim(inits$sel_coff)[1])
  inits$sel_coff[9:11,2,1:8] <- rep(log_selcoffs2, each = 3) # Males

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)

  # Check selectivity is 0 for bins 1-2
  output <- c(ss_run$quantities$sel_at_age[,1,1:2,])
  testthat::expect_equal(output, rep(0, length(output)))
})
