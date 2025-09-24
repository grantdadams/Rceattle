library(Rceattle)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060
BS2017SS$fleet_control$proj_F_prop <-rep(1,7)


# Harvest projections ----
# Single-species with fixed M
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)


# -- NPFMC Tier 3
ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                  inits = ss_run$estimated_params, # Initial parameters from ss_run
                                  estimateMode = 2, # Run projection only
                                  HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                  Ftarget = 0.4, # F40%
                                                  Flimit = 0.35, # F35%
                                                  Plimit = 0.2, # No fishing when SB<SB20
                                                  Alpha = 0.05),
                                  msmMode = 0, # Single species mode
                                  verbose = 1)

# -- NPFMC Tier 3 (dynamic B0)
ss_run_dynamicTier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 2, # Run projection only
                                         HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         Ftarget = 0.4, # F40%
                                                         Flimit = 0.35, # F35%
                                                         Plimit = 0.2, # No fishing when SB<SB20
                                                         Alpha = 0.05),
                                         msmMode = 0, # Single species mode
                                         verbose = 1)


# Environmental projections ----
# R = exp(logMu + B * X + e_y)
# logR = logMu + B * X + e_y
# (logR - logMu - e_y)/B = X

# - Create random env data
set.seed(123)
nyrs <- length(BS2017SS$styr:BS2017SS$endyr)
nyrsp <- length(BS2017SS$styr:BS2017SS$projyr)
projyrs <- (BS2017SS$endyr + 1):BS2017SS$projyr
e_y <- rnorm(nyrs)
beta <- c(10, 0.2, -.3)
env_data <- data.frame(Year = BS2017SS$styr:BS2017SS$projyr)
for(i in 1:3){
  env_data[paste0("Ind",i)]= c((log(ss_run$quantities$R[i,1:nyrs]) - log(mean(ss_run$quantities$R[i,1:nyrs]) - e_y))/beta[i], rnorm(length(projyrs))) + rnorm(nyrsp, 0, 0.2)
}

plot(env_data$Ind1, ss_run$quantities$R[1,])

BS2017SS$env_data <- BS2017SS$env_data %>% full_join(env_data)


# Three species model (indX impacts speciesX recruitment)
sem_3spp = "
  # link, lag, param_name, start_value
  recdevs1 <-> recdevs1, 0, sigmaR1, 1
  recdevs2 <-> recdevs2, 0, sigmaR2, 1
  recdevs3 <-> recdevs3, 0, sigmaR3, 1
  Ind1 -> recdevs1, 0, ind1_to_R1, 0
  Ind2 -> recdevs2, 0, ind2_to_R2, 0
  Ind3 -> recdevs3, 0, ind3_to_R3, 0
"

proj_model <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = TRUE, # Random recruitment
                            dsem = build_DSEM(
                              sem = sem_3spp,
                              family = "normal",
                              estimate_projection = TRUE,
                            ),
                            msmMode = 0, # Single species mode
                            phase = FALSE,
                            verbose = 1)

check <- data.frame(Parname = names(proj_model$opt$opt$par), Par = proj_model$opt$opt$par, Gr = proj_model$obj$gr())

summ <- summary(proj_model$dsem)
summ$Estimate <- proj_model$estimated_params$beta_z
summ


plot_ssb(list(ss_run, proj_model), model_names = c("Default", "DSEM"), incl_proj = TRUE)
plot_recruitment(list(ss_run, proj_model), model_names = c("Default", "DSEM"), incl_proj = TRUE)
proj_model$quantities$rec_dev

