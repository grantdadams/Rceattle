library(Rceattle)

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
data("EBS_ss_run")
BS2017SS$projyr <- 2100

# Write data to excel
Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")


################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
mydata$fleet_control$proj_F_prop <-rep(1,7)
params <- build_params(mydata)
inits <- EBS_ss_run$estimated_params
inits$log_gam_a <- params$log_gam_a
inits$log_gam_b <- params$log_gam_b
inits$log_phi <- params$log_phi
inits$logH_1 <- params$logH_1
inits$logH_1a <- params$logH_1a
inits$logH_1b <- params$logH_1b
inits$logH_2 <- params$logH_2
inits$logH_3 <- params$logH_3
inits$H_4 <- params$H_4
inits$ln_FSPR <- NULL
inits$ln_Ftarget <- params$ln_Ftarget
inits$ln_Flimit <- params$ln_Flimit

ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default",
                            verbose = 1)
length(ss_run$opt$par)
length(EBS_ss_run$opt$par)

plot_biomass(list(ss_run, EBS_ss_run), model_names = c(1,2), incl_proj = TRUE)
plot_catch(list(ss_run, EBS_ss_run), model_names = c(1,2), incl_proj = TRUE)
plot_logindex(list(ss_run, EBS_ss_run), model_names = c(1,2), incl_proj = TRUE)

sum(ss_run$quantities$sel - EBS_ss_run$quantities$sel)
sum(ss_run$quantities$jnll_comp[1:13,] != EBS_ss_run$quantities$jnll_comp[1:13,])
round(ss_run$quantities$jnll_comp[1:13,] - EBS_ss_run$quantities$jnll_comp[1:13,],3)
ss_run$quantities$jnll_comp[,1:3]
EBS_ss_run$quantities$jnll_comp[1:18,1:3]


plot_biomass(list(ss_run, EBS_ss_run), model_names = c(1,2), incl_proj = TRUE)
plot_catch(list(ss_run, EBS_ss_run))
plot_selectivity(ss_run);
plot_selectivity(EBS_ss_M_run)
ss_run$identified$BadParams[which(ss_run$identified$BadParams$Param_check == "Bad"),]
ss_run$identified$param_list$rec_dev
