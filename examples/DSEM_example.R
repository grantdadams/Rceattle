# Code to run DSEM-linked assessment models
# NOTE: we have to phase the models to get to a better starting point
library(Rceattle)


# Model 1 ----
# Single-species GOA pollock, IID recruitment
data("GOApollock") # Single-species data. ?BS2017SS for more information on the data
GOApollock$projyr <- 2020  # Shorten proj year for quicker estimation
GOApollock$env_data <- GOApollock$env_data %>%
  dplyr::mutate(ScaledBT =scale(BTempC))

model1 <- Rceattle::fit_mod(data_list = GOApollock,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = TRUE, # Random recruitment
                            msmMode = 0, # Single species mode
                            phase = FALSE,
                            verbose = 1)
model1$data_list$dsem_settings$sem


# Model 2 ----
# IID sem, if start value is not provided, DSEM terms will not be estimated
sem_iid = "
  # link, lag, param_name, start_value

  #internal data relation
  ScaledBT -> ScaledBT, 1, AR_BT, 0
  ScaledBT -> recdevs1, 1, BT_to_R, 0
  recdevs1 <-> recdevs1, 0, sigmaR1, 1
"
model2 <- Rceattle::fit_mod(data_list = GOApollock,
                            inits = NULL, # Initial parameters from model 1
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            dsem = build_DSEM(
                              sem = sem_iid,
                              family = "normal"
                            ),
                            random_rec = TRUE, # Random recruitment
                            msmMode = 0, # Single species mode
                            phase = FALSE,
                            verbose = 1)


# * Plot ----
mod_list <- list(model1, model2)
mod_names <- c("Single-species IID", "Single-species Env")
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names)


# Model 3 ----
# Three species IID model
data("BS2017SS") # Single-species data for EBS model
BS2017SS$projyr <- 2020  # Shorten proj year for quicker estimation
model3 <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = TRUE, # Random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)
print(model3$data_list$dsem_settings$sem) # Default SEM
model3$estimated_params$beta_z # sigmaR

# Model 4 ----
# Three species IID model (cod and ATF recruitment impact pollock recruitment)
sem_3spp = "
  # link, lag, param_name, start_value
  recdevs1 <-> recdevs1, 0, sigmaR1, 1
  recdevs2 <-> recdevs2, 0, sigmaR2, 1
  recdevs3 <-> recdevs3, 0, sigmaR3, 1
  recdevs2 -> recdevs1, 0, R2_to_R1, 0
  recdevs3 -> recdevs1, 0, R3_to_R1, 0
"

model4 <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = TRUE, # Random recruitment
                            dsem = build_DSEM(
                              sem = sem_3spp,
                              family = "normal"
                            ),
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

# * Plot ----
mod_list <- list(model3, model4)
mod_names <- c("Three-species IID", "Three-species linked")
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names)

