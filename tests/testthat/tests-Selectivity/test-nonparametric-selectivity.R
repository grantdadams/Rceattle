
testthat::test_that("Test age-based non-parametric selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Adjust data
  GOA2018SS$fleet_control$Selectivity <- 2
  GOA2018SS$fleet_control$N_sel_bins <- 8
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_curve_pen1 <- 5
  GOA2018SS$fleet_control$Sel_curve_pen2 <- 10
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1

  inits <- suppressMessages(build_params(GOA2018SS))


  # ADMB code
  # dvar_vector sel_coffs_tmp(1,n_sel_bins_fsh(k));
  # for (i=styr;i<=endyr;i++)
  # {
  #   if (i==yrs_sel_ch_fsh(k,isel_ch_tmp))
  #   {
  #     sel_coffs_tmp.initialize();
  #     sel_coffs_tmp = log_selcoffs_fsh(k,isel_ch_tmp);
  #     avgsel_fsh(k,isel_ch_tmp)              = log(mean(mfexp(sel_coffs_tmp)));
  #   }
  #   log_sel_fsh(k,i)(1,n_sel_bins_fsh(k))        = sel_coffs_tmp;
  #   log_sel_fsh(k,i)(n_sel_bins_fsh(k),nages)    = log_sel_fsh(k,i,n_sel_bins_fsh(k));
  #   log_sel_fsh(k,i)                                  -= log(mean(mfexp(log_sel_fsh(k,i) )));
  # }

  # Specify selectivity
  n_sel_bins <- 8
  log_selcoffs <- rnorm(n_sel_bins)
  log_selcoffs2 <- rnorm(n_sel_bins)
  avgsel_fsh <- log(mean(exp(log_selcoffs)))

  # - Pk
  log_sel_fsh <- c(log_selcoffs, rep(log_selcoffs[n_sel_bins], GOA2018SS$nages[1]-n_sel_bins))
  log_sel_fsh <- log_sel_fsh - log(mean(exp(log_sel_fsh)))

  # - ATF
  log_sel_fsh2 <- c(log_selcoffs, rep(log_selcoffs[n_sel_bins], GOA2018SS$nages[2]-n_sel_bins))
  log_sel_fsh2 <- log_sel_fsh2 - log(mean(exp(log_sel_fsh2)))

  log_sel_fsh2m <- c(log_selcoffs2, rep(log_selcoffs2[n_sel_bins], GOA2018SS$nages[2]-n_sel_bins))
  log_sel_fsh2m <- log_sel_fsh2m - log(mean(exp(log_sel_fsh2m)))

  # - Cod
  log_sel_fsh3 <- c(log_selcoffs, rep(log_selcoffs[n_sel_bins], GOA2018SS$nages[3]-n_sel_bins))
  log_sel_fsh3 <- log_sel_fsh3 - log(mean(exp(log_sel_fsh3)))

  # Set params
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

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh), tolerance = 0.0001))

  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh2), tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh2m), tolerance = 0.0001))

  # - Cod
  apply(ss_run$quantities$sel_at_age[12,1,1:12,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh3), tolerance = 0.0001))

  # Check penalties
  # - Curvature penalty
  pen1 <- sum(10*(diff(diff(log_sel_fsh))^2))

  # - Descending penalty
  difftmp = -diff(log_sel_fsh) # Decreasing will be positive
  difftmp = (abs(difftmp) + difftmp)/2
  pen2 <- sum(difftmp^2) * 5

  # - Avgsel
  pen3 <- 2 * log(mean(exp(log_selcoffs)))^2

  testthat::expect_equal(pen1+pen2+pen3, as.numeric(ss_run$quantities$jnll_comp[5,8]), tolerance = 0.0001)
})

testthat::test_that("2DAR1 selectivity map and likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  GOA2018SS$fleet_control$Selectivity <- 6 # Age-based 2DAR1
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$N_sel_bins <- 8
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = NULL, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = TRUE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      phase = TRUE,
                      verbose = 0)

  # Hyper parameters (on except for not estimated fleet)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_ln_sd), c(1:6, NA, 8:16))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,1]), c(1:6, NA, 8:16))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,2]), c(1:6, NA, 8:16) + 16)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,3]), as.numeric(rep(NA, 16)))

  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_dev_ln_sd), rep(0, 16))
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_curve_pen[,1]), rep(0, 16))
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_curve_pen[,2]), rep(0, 16))

  # Fixed effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,,1:8]), rep(0, 2 * 8 * 16))
  # - Females (sex combined)
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[-7,1,1:8])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[7,1,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[-7,1,1:8])), 8 * 15)

  # - Males
  flt2sex <- GOA2018SS$nsex[GOA2018SS$fleet_control$Species] == 2
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[!flt2sex,2,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])), 8 * sum(flt2sex))

  # Random effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,1,1:8]), rep(0, 8 * 16))

  # - Females (sex combined)
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[-7,1,1:8,])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[7,1,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[-7,1,1:8,])), 8 * 15 * nyrs)


  # - Males
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[!flt2sex,2,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])), 8 * sum(flt2sex) * nyrs)

  # TMB object
  testthat::expect_equal(length(unique(ss_run$obj$env$random)),  8 * 18 * nyrs)

  # JNLL
  testthat::expect_all_true(ss_run$quantities$jnll_comp[6,-7] != 0)
  testthat::expect_all_true(ss_run$quantities$jnll_comp[6,7] == 0)
  testthat::expect_all_true(!is.na(ss_run$quantities$jnll_comp[6,]))
})

testthat::test_that("3DAR1 selectivity map and likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  GOA2018SS$fleet_control$Selectivity <- 7 # Age-based 3DAR1
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel_sd_prior <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$N_sel_bins <- 8
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              random_sel = TRUE, # Turn on laplace for sel devs
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 0)

  # Hyper parameters (on except for not estimated fleet)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_ln_sd), c(1:6, NA, 8:16))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,1]), c(1:6, NA, 8:16))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,2]), c(1:6, NA, 8:16) + 16)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,3]), c(1:6, NA, 8:16) + 2*16)

  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_dev_ln_sd), rep(0, 16))
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_curve_pen[,1]), rep(0, 16))
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_curve_pen[,2]), rep(0, 16))

  # Fixed effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,,1:8]), rep(0, 2 * 8 * 16))
  # - Females (sex combined)
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[-7,1,1:8])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[7,1,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[-7,1,1:8])), 8 * 15)

  # - Males
  flt2sex <- GOA2018SS$nsex[GOA2018SS$fleet_control$Species] == 2
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[!flt2sex,2,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])), 8 * sum(flt2sex))

  # Random effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,1,1:8]), rep(0, 8 * 16))

  # - Females (sex combined)
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[-7,1,1:8,])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[7,1,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[-7,1,1:8,])), 8 * 15 * nyrs)


  # - Males
  testthat::expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])))
  testthat::expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[!flt2sex,2,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])), 8 * sum(flt2sex) * nyrs)

  # TMB object
  testthat::expect_equal(length(unique(ss_run$obj$env$random)),  8 * 18 * nyrs)

  # JNLL
  testthat::expect_all_true(ss_run$quantities$jnll_comp[6,-7] != 0)
  testthat::expect_all_true(ss_run$quantities$jnll_comp[6,7] == 0)
  testthat::expect_all_true(!is.na(ss_run$quantities$jnll_comp[6,]))
})
