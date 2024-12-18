data("GOA2018SS") # ?BS2017SS for more information on the data
GOA2018SS$projyr <- 2030

ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                            file = NULL,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1,
                            phase = "default")

ss_runre <- Rceattle::fit_mod(data_list = GOA2018SS,
                            file = NULL,
                            inits = ss_run$estimated_params, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = TRUE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1,
                            phase = NULL)

ss_run_se <- Rceattle::fit_mod(data_list = GOA2018SS,
                            file = NULL,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1,
                            random_sel = TRUE,
                            phase = "default")




check <- list()

for(i in 1:length(ss_run$map$mapList)){
  map <- ss_run$map
  map$mapFactor[[i]] <- factor(map$mapList[[i]] * NA)
  check[[i]] <- try(Rceattle::fit_mod(data_list = BS2017SS,
                                      file = NULL,
                                      map = map,
                                      inits = ss_run$estimated_params, # Initial parameters = 0
                                      estimateMode = 0, # Estimate
                                      random_rec = TRUE, # No random recruitment
                                      msmMode = 0, # Single species mode
                                      verbose = 1,
                                      phase = NULL), silent = TRUE)
}


