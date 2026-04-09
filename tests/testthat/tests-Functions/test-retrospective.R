testthat::test_that("Test retrospective", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)


  # Data
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$projyr <- 2020
  BS2017SS$fleet_control$proj_F_prop <-rep(1,7)


  # Single-species with fixed M
  ss_run <- fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 1, # Estimate hindcast only
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

  # Retro
  ret <-retrospective(ss_run, peels = 5)

  # By hand
  retro_list <- list()
  for(y in 1:5){
    dat_tmp <- BS2017SS
    dat_tmp$endyr <- dat_tmp$endyr - y
    retro_list[[y]] <- Rceattle::fit_mod(data_list = dat_tmp,
                                         inits = NULL, # Initial parameters = 0
                                         file = NULL, # Don't save
                                         estimateMode = 1, # Estimate hindcast only
                                         random_rec = FALSE, # No random recruitment
                                         msmMode = 0, # Single species mode
                                         phase = TRUE,
                                         loopnum = 5,
                                         verbose = 1)
  }
  retro_list <- rev(retro_list)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Tests ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  yrs <- ss_run$data_list$styr:ss_run$data_list$endyr
  nyrs <- length(yrs)

  testthat::expect_equal(6, length(ret$Rceattle_list))
  testthat::expect_equal(2012:2017, sapply(ret$Rceattle_list, function(x) x$data_list$endyr_peel))

  #testthat::expect_equal(lapply(retro_list, function(x) x$quantities$biomass[,1:nyrs]), lapply(ret$Rceattle_list[1:5], function(x) x$quantities$biomass[,1:nyrs]), tolerance = 0.01)
})
