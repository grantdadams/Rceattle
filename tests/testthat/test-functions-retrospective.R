# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test retrospective", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)


  # Data
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$projyr <- 2020
  BS2017SS$fleet_control$Proj_F_proportion <-rep(1,7)


  # Single-species with fixed M
  ss_run <- fit_mod(data_list = BS2017SS,
                    inits = NULL, # Initial parameters = 0
                    file = NULL, # Don't save
                    estimateMode = 1, # Estimate hindcast only
                    random_rec = FALSE, # No random recruitment
                    msmMode = 0, # Single species mode
                    fit_control = fit_control(getsd = FALSE, 
                      phase = TRUE,
                      verbose = 1))

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
                                         fit_control = fit_control(getsd = FALSE, 
                                           phase = TRUE,
                                           loopnum = 5,
                                           verbose = 1))
  }
  retro_list <- rev(retro_list)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Tests ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  yrs <- ss_run$data_list$styr:ss_run$data_list$endyr
  nyrs <- length(yrs)

  testthat::expect_equal(6, length(ret$Rceattle_list))
  testthat::expect_equal(2012:2017, as.numeric(sapply(ret$Rceattle_list, function(x) x$data_list$endyr_peel)))

  for(i in 1:5){
    # Matches
    yrs <- retro_list[[i]]$data_list$styr:retro_list[[i]]$data_list$endyr
    nyrs <- length(yrs)
    testthat::expect_equal(retro_list[[i]]$quantities$biomass[,1:nyrs], ret$Rceattle_list[[i]]$quantities$biomass[,1:nyrs], tolerance = 0.01)

    # Endyr of peel
    testthat::expect_equal(ret$Rceattle_list[[i]]$data_list$endyr_peel, 2017-rev(1:5)[i])

    # Data removed
    # * Index
    index_tmp <- ret$Rceattle_list[[i]]$data_list$index_data %>% filter(Year > 2017-rev(1:5)[i])
    testthat::expect_equal(nrow(index_tmp), 0)

    # * Comp
    comp_tmp <- ret$Rceattle_list[[i]]$data_list$comp_data %>% filter(Year > 2017-rev(1:5)[i])
    testthat::expect_equal(nrow(comp_tmp), 0)

    # * CAAL
    caal_tmp <- ret$Rceattle_list[[i]]$data_list$caal_data %>% filter(Year > 2017-rev(1:5)[i])
    testthat::expect_equal(nrow(caal_tmp), 0)
  }
})
