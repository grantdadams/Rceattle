## Load packages ----
# install.packages("pacman")
pacman::p_load(forecast)

## 1) Simulate data ----
# * Define population and observation model ----
# Mortality and population specifications
M <- 0.2
nyrs <- 50
yrs <- 1:nyrs
nages <- 10
ages <- 1:nages

# - Weight at age
winf <- 10
vb_k <- 0.25
wt <- winf * (1-exp(-vb_k * ages))

# - Maturity
mat <- as.numeric(ages>4) # Knife-edge maturity above age 4

# - Recruitment
Rhat <- 1000
Rdev <- arima.sim(list(order=c(1,0,0), ar=.5), n = nyrs, rand.gen = rlnorm) # Follows AR1 process
R <- Rhat * Rdev
plot(y = R, x = yrs, type = "l")

# - Fishing mortality and selectivity
# -- here I am assuming the seperability where there is an increasing and decreasing fishing mortality
# -- and logistic fishery selectivity
fsh_mort <- c(seq(0, 0.3, length.out = nyrs-15), seq(0.3, 0.1, length.out = 15)) # Increasing and decreasing F
fsh_slope <- 2
fsh_asymp <- 4
fsh_sel <- 1/(1+exp(-fsh_slope * (ages-fsh_asymp)))
plot(x = ages, y = fsh_sel, type = "l")

# - Survey selectivity and catchability
srv_q <- 0.4
srv_slope <- 2
srv_asymp <- 3
srv_sel <- 1/(1+exp(-srv_slope * (ages-srv_asymp))) # Logistic selectivity that targets young fish
plot(x = ages, y = srv_sel, type = "l")


# * Fill population and observation model ----
n_at_age <- matrix(0, nrow = nages, ncol = nyrs + 1)
colnames(n_at_age) <- 0:nyrs; rownames(n_at_age) <- ages
srv_comp <- fsh_comp <- matrix(0, nrow = nages, ncol = nyrs)
biomass <- index_data <- ssb <- catch <- c()

# - Calculate initial abundance and recruitment
n_at_age[1,2:(nyrs+1)] <- R # Age-1 is recruits
n_at_age[1,1] <- Rhat # Age-1, Yr-0 is Rhat
n_at_age[2:nages,1] <- Rhat * exp(-M * 1:(nages-1))
n_at_age[nages,1] <- n_at_age[nages,1]/(1-exp(-M*nages)) # Max-age, solve for geometric series for equilibrium

# - Loop through years and ages
for(y in 2:nyrs){ # Year one starts at equilibrium
  for(a in 2:nages){
    n_at_age[a,y] <- n_at_age[a-1,y-1] * exp(-M - fsh_sel[a-1] * fsh_mort[y-1])
  }

  # - Calculate biomass
  biomass[y] <- sum(n_at_age[,y] * wt)
  ssb[y] <- sum(n_at_age[,y] * wt * mat * 0.5) # Assuming 50% females


  # - Calculate fishery and survey age composition
  srv_comp[,y] <- n_at_age[,y] * srv_sel # Assuming happens at month-0
  fsh_comp[,y] <- n_at_age[,y] * (fsh_sel * fsh_mort[y])/ (M + fsh_sel * fsh_mort[y]) * (1-exp(-M - fsh_sel * fsh_mort[y])) # Baranov catch equation
  # - Baranov catch equation is the integral of dN/dT -(F + M) times the ratio of F to (M+F)
  # - The 1 - exp(-M+F) is calculating the number of fish lost due to F + M: (N-N*exp(-M+F))

  # - Calculate catch and survey biomass
  index_data[y] <- sum(srv_comp[,y] * wt * srv_q)
  catch[y] <- sum(fsh_comp[,y] * wt)

  # - Normalize comp data
  srv_comp[,y] <- srv_comp[,y]/sum(srv_comp[,y])
  fsh_comp[,y] <- fsh_comp[,y]/sum(fsh_comp[,y])
}


# * Add observation error ----
catch_logsd <- 0.05
srv_logsd <- 0.2
fsh_compn <- 1000
srv_compn <- 1000

# - Simulate observed catch and survey data from lognormal distribution
srv_obs <- rlnorm(n = nyrs, meanlog = log(index_data), sd = srv_logsd)
catch_obs <- rlnorm(n = nyrs, meanlog = log(catch), sd = catch_logsd)
catch_obs[is.na(catch_obs)] <- 0

# - Simulate age composition data from multinomial with assumed sample size of 200 for survey and 100 for fishery
# -- Size is the multinomial sample size (See Hamel and Stewart 2014 for calculation with observed data)
srv_comp_obs <- apply(srv_comp, 1, function(x) rmultinom(n = 1, size = srv_compn, prob = x))
srv_comp_obs <- apply(srv_comp_obs, 1, function(x) x/sum(x)) # Normalize

fsh_comp_obs <- apply(fsh_comp, 1, function(x) rmultinom(n = 1, size = fsh_compn, prob = x))
fsh_comp_obs <- apply(fsh_comp_obs, 1, function(x) x/sum(x))

# * Plot ----
plot(x = 1:nyrs, y = biomass, ylab = "Biomass", type = "l", col = 2)

# - Survey
plot(x = 1:nyrs, y = srv_obs, ylab = "Srv data", type = "l", col = 2)
lines(x = 1:nyrs, y = index_data, lty = 2, col = 2)

# - Catch
plot(x = 1:nyrs, y = catch_obs, ylab = "Catch", type = "l", col = 1)
lines(x = 1:nyrs, y = catch, lty = 2, col = 1)

## 2) Convert to Rceattle data ----
library(Rceattle)
data("GOApollock")
sim_dat <- GOApollock

# * Control ----
sim_dat$nspp <- 1
sim_dat$styr <- 1
sim_dat$endyr <- nyrs
sim_dat$projyr <- nyrs + 1
sim_dat$nsex <- 1
sim_dat$spawn_month <- 0
sim_dat$R_sexr <- 0.5
sim_dat$nages <- nages
sim_dat$minage <- 1
sim_dat$nlengths <- 10
sim_dat$pop_wt_index <- 1
sim_dat$ssb_wt_index <- 1
sim_dat$pop_age_transition_index <- 1
sim_dat$sigma_rec_prior <- 1.24
sim_dat$other_food <- 1e6
sim_dat$estDynamics <- 0
sim_dat$est_sex_ratio <- 0
sim_dat$sex_ratio_sigma <- 1

# * Fleet control ----
sim_dat$fleet_control <-
  data.frame(Fleet_name = c("Survey", "Fishery"),
             Fleet_code = 1:2,           # 1) Temporary survey index
             Fleet_type = 2:1,           # 2) Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
             Species = 1,              # 3) Species
             Selectivity_index = 1:2,    # 4) Survey selectivity index
             Selectivity = 1,          # 5) Selectivity type
             N_sel_bins = NA,             # 6) Non-parametric selectivity ages
             Time_varying_sel = 0,     # 7) Time-varying selectivity type.
             Time_varying_sel_sd_prior = 0,
             Bin_first_selected = 1,   # 8) First age selected
             Sel_norm_bin1 = NA,       # 9b) Age of max selectivity (used for normalization). If NA, does not normalize
             Sel_norm_bin2 = NA,       # 9a) upper age of max selectivity (used for normalization). If NA, does not normalize
             Comp_loglike = 0,         # 10) Index indicating wether to do dirichlet multinomial for a multinomial)
             Weight1_Numbers2 = 1,     # 11) Survey units
             Weight_index = 1,         # 12) Dim1 of weight (what weight-at-age data set)
             Age_transition_index = 1, # 13) Dim3 of age transition matrix (what ALK to use)
             Q_index = c(1, NA),              # 14) Index of survey q
             Catchability = c(1,0),           # 15) Parametric form of q
             Q_prior = 1,
             Q_sd_prior = 0.1,
             Time_varying_q = 0,       # 16) Time varying q type
             Time_varying_q_sd_prior = 0,
             Estimate_index_sd = 0,    # 17) Wether to estimate standard deviation of survey time series
             Index_sd_prior = NA,
             Estimate_catch_sd = 0,    # 18) Wether to estimate standard deviation of fishery time series)
             Catch_sd_prior = NA,
             proj_F_prop  = 0,
             Comp_weights = 1)

# * Survey data ----
# - Starts at year-2
sim_dat$index_data <- data.frame(
  Fleet_name = "Survey",
  Fleet_code = 1,
  Species = 1,
  Year = 2:sim_dat$endyr,
  Month = 0,
  Selectivity_block = 1,
  Q_block = 1,
  Observation = srv_obs[2:sim_dat$endyr],
  Log_sd = srv_logsd
)

# * Fishery data ----
# - Starts at year-2
sim_dat$catch_data <- data.frame(
  Fleet_name = "Fishery",
  Fleet_code = 2,
  Species = 1,
  Year = sim_dat$styr:sim_dat$endyr,
  Month = 0,
  Selectivity_block = 1,
  Catch = catch_obs,
  Log_sd = catch_logsd
)


# * Age comp ----
# - Survey
comp_data_srv <- data.frame(
  Fleet_name = "Survey",
  Fleet_code = 1,
  Species = 1,
  Sex = 0,
  Age0_Length1 = 0,
  Year = 2:sim_dat$endyr,
  Month = 0,
  Sample_size = 200) # Set as above

comp_data_srv <- cbind(comp_data_srv, t(srv_comp_obs)[2:sim_dat$endyr, ])

# - Fishery
comp_data_fsh <- data.frame(
  Fleet_name = "Fishery",
  Fleet_code = 2,
  Species = 1,
  Sex = 0,
  Age0_Length1 = 0,
  Year = 2:sim_dat$endyr,
  Month = 0,
  Sample_size = 100) # Set as above

comp_data_fsh <- cbind(comp_data_fsh, t(fsh_comp_obs)[2:sim_dat$endyr, ])

# - Combine
comp_data <- rbind(comp_data_srv, comp_data_fsh)
colnames(comp_data) <- colnames(sim_dat$comp_data)[1:(8+nages)]
sim_dat$comp_data <- comp_data

# * Emp sel and nbyage ----
sim_dat$emp_sel <- sim_dat$emp_sel[0,]
sim_dat$NByageFixed <- sim_dat$NByageFixed[0,]

# * Age transition matrix (age to length) ----
sim_dat$age_trans_matrix <- sim_dat$age_trans_matrix[0,1:15]
age_trans_matrix <- data.frame(
  Age_transition_name = "Species1",
  Age_transition_index = 1,
  Species = 1,
  Sex = 0,
  Age = 1:sim_dat$nages
)
age_trans_matrix <- cbind(age_trans_matrix, diag(nages))
colnames(age_trans_matrix) <- colnames(sim_dat$age_trans_matrix)
sim_dat$age_trans_matrix  <- age_trans_matrix


# * Age error matrix ----
sim_dat$age_error <- sim_dat$age_error[0,1:12]
age_error <- data.frame(
  Species = 1,
  True_age = 1:sim_dat$nages
)
age_error <- cbind(age_error, diag(nages))
colnames(age_error) <- colnames(sim_dat$age_error)
sim_dat$age_error  <- age_error


# * Weight-at-age ----
sim_dat$weight <- sim_dat$weight[0,1:(nages+5)]
sim_dat$weight[1,] <-c("Species1", 1, 1, 0 ,0, wt)


# * Maturity ----
sim_dat$maturity <- sim_dat$maturity[0,1:(nages+1)]
sim_dat$maturity[1,] <-c(1, mat)


# * Sex-ratio ----
sim_dat$sex_ratio <- sim_dat$sex_ratio[0,1:(nages+1)]
sim_dat$sex_ratio[1,] <-c(1, rep(0.5, nages))


# * Mortality ----
sim_dat$M1_base <- sim_dat$M1_base[0,1:(nages+2)]
sim_dat$M1_base[1,] <-c(1, 0, rep(M, nages))

# * Envdata ----
sim_dat$env_data <- data.frame(
  Year = 1:nyrs,
  BTempC = 5
)

# * ration_data and diet ----
sim_dat$ration_data <- sim_dat$ration_data[0,]
sim_dat$diet_data <- sim_dat$diet_data[0,]


## 3) Fit CEATTLE ----
sim_model <- Rceattle::fit_mod(
  data_list = sim_dat,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = FALSE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = TRUE) # Use default phasing

plot(x = 1:nyrs, y = sim_model$quantities$biomass[1,1:nyrs], ylab = "Biomass", type = "l", col = 2)
lines(x = 1:nyrs, y = biomass)

