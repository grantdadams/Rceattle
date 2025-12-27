
test_that("non-parametric selectivity not normalized", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Adjust data
  GOA2018SS$fleet_control$Selectivity <- 2
  GOA2018SS$fleet_control$N_sel_bins <- 8
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_curve_pen1 <- 5
  GOA2018SS$fleet_control$Sel_curve_pen2 <- 10
  GOA2018SS$fleet_control$Sel_norm_bin1 <- NA
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Sel_sd_prior <- 1

  inits <- build_params(GOA2018SS)


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
                              verbose = 1)

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh), tolerance = 0.0001))

  # - ATF
  apply(ss_run$quantities$sel[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh2), tolerance = 0.0001))
  apply(ss_run$quantities$sel[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh2m), tolerance = 0.0001))

  # - Cod
  apply(ss_run$quantities$sel[12,1,1:12,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh3), tolerance = 0.0001))

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
