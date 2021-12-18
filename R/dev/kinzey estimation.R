library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$est_M1 <- c(0,0,0)


################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
BS2017SS$fleet_control$proj_F_prop <- c(rep(1,3), rep(0,4))
BS2017SS$projyr = 2100
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = "EBS_ss_run", # Don't save
                            debug = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1)
data("EBS_ss_run")
plot_biomass(list(EBS_ss_run, ss_run))


data("BS2017MS")
data("EBS_ms_run")
load("~/GitHub/Rceattle/EBS_ms_run_1iter.RData")
BS2017MS$projyr = 2100
inits <- mod_objects$estimated_params
inits$log_gam_a <- ss_run$estimated_params$log_gam_a
inits$log_gam_b <- ss_run$estimated_params$log_gam_b
inits$log_phi <- ss_run$estimated_params$log_phi
inits$logH_1 <- ss_run$estimated_params$logH_1
inits$logH_1a <- ss_run$estimated_params$logH_1a
inits$logH_1b <- ss_run$estimated_params$logH_1b
inits$logH_2 <- ss_run$estimated_params$logH_2
inits$logH_3 <- ss_run$estimated_params$logH_3
inits$H_4 <- ss_run$estimated_params$H_4

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            debug = 0, # Estimate
                            niter = 3, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)
# Check starting nll
for(i in 1:4){
  ms_run_check <- Rceattle::fit_mod(data_list = BS2017MS,
                                    inits = ms_run$estimated_params, # Initial parameters from single species ests
                                    file = NULL, # Don't save
                                    debug = 1, # Estimate
                                    niter = 3, # 10 iterations around population and predation dynamics
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = 1, # MSVPA based
                                    suitMode = i, # empirical suitability
                                    verbose = 1)
  plot_biomass(ms_run_check)
  print(i)
  print(ms_run_check$quantities$jnll_comp)
  save(i,file = paste0("Suit",i,".RData"))
}

for(i in 3:9){ # Check 8
  ms_run_check <- Rceattle::fit_mod(data_list = BS2017MS,
                                    inits = ms_run$estimated_params, # Initial parameters from single species ests
                                    file = NULL, # Don't save
                                    debug = 1, # Estimate
                                    niter = 3, # 10 iterations around population and predation dynamics
                                    random_rec = FALSE, # No random recruitment
                                    msmMode = i, # MSVPA based
                                    suitMode = 1, # empirical suitability
                                    verbose = 1)
  plot_biomass(ms_run_check)
  print(i)
  print(ms_run_check$quantities$jnll_comp)
  save(i,file = paste0("MsmMode",i,".RData"))
}
sum(ms_run_check$quantities$ration_hat)
sum(ms_run_check$quantities$ration_hat_ave)
for(i in 1:42){
  print(sum((ms_run_check$quantities$ration_hat[1,1,1:12,i] - ms_run_check$quantities$ration_hat_ave[1,1,1:12])^2))
}
ms_run_check$quantities$B_eaten_as_pred[1,1,1:12,1:42]/ms_run_check$quantities$AvgN[1,1,1:12,1:42]


# Try estimating
BS2017MS$est_M1 <- rep(0,3)
ms_suit <- list()
for(i in 1:4){
  ms_suit[[i]] <- list()
  for(j in 3:9){
    ms_suit[[i]][[j]] <- try(Rceattle::fit_mod(data_list = BS2017MS,
                                          inits = ss_run$estimated_params, # Initial parameters from single species ests
                                          file = NULL, # Don't save
                                          debug = 0, # Estimate
                                          niter = 3, # 10 iterations around population and predation dynamics
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = i, # MSVPA based
                                          suitMode = j, # empirical suitability
                                          verbose = 1))
  }
}

plot_biomass(c(list(ss_run),ms_suit[[1]][5:9], ms_suit[[2]][5:9]))
