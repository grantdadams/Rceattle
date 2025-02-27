# TODO: fleet_code to fleet_index
# Get rid of sex ratio in control
# call "wt" "weight"
# call "pmature" "maturity"

# Simulate Data -----------------------------------------------------------
# From Matt Cheng
sim_pop_model <- function(years,
                          ages,
                          WAA,
                          MatAA,
                          mean_Rec,
                          sigma_R,
                          sigma_catch,
                          sigma_srv,
                          fish_ISS,
                          srv_ISS,
                          M,
                          fish_sel,
                          srv_sel,
                          Fmort,
                          srv_q) {

  n_yrs <- length(years)
  n_ages <- length(ages)

  # Initialize arrays
  NAA <- array(0, dim = c(n_yrs + 1, n_ages))  # Numbers at age
  ZAA <- array(0, dim = c(n_yrs, n_ages))      # Total mortality
  CAA <- array(0, dim = c(n_yrs, n_ages))      # Catch at age
  FAA <- array(0, dim = c(n_yrs, n_ages))      # Fishing mortality at age

  # Population metrics
  SSB <- numeric(n_yrs)
  Total_Biom <- numeric(n_yrs)
  Catch <- numeric(n_yrs)

  # Generate recruitment deviations
  rec_devs <- rnorm(n_yrs, 0, sigma_R)
  init_devs <- rnorm(n_ages - 2, 0, sigma_R)

  # Initialize population
  init_age_idx <- 1:(n_ages - 2)
  NAA[1, init_age_idx + 1] <- mean_Rec * exp(init_devs - (init_age_idx * M))
  NAA[1, n_ages] <- mean_Rec * exp(-(n_ages - 1) * M) / (1 - exp(-M))

  # Project population forward
  for(y in 1:n_yrs) {
    # New recruits
    NAA[y, 1] <- mean_Rec * exp(rec_devs[y])

    # Calculate mortality
    FAA[y,] <- Fmort[y] * fish_sel
    ZAA[y,] <- FAA[y,] + M

    # Calculate catch
    CAA[y,] <- FAA[y,] / ZAA[y,] * NAA[y,] * (1 - exp(-ZAA[y,]))

    # Project survivors
    if(y < n_yrs) {
      for(a in 1:(n_ages-1)) {
        NAA[y+1, a+1] <- NAA[y, a] * exp(-ZAA[y, a])
      }
      # Plus group
      NAA[y+1, n_ages] <- NAA[y+1, n_ages] + NAA[y, n_ages] * exp(-ZAA[y, n_ages])
    }

    # Calculate annual metrics
    Total_Biom[y] <- sum(NAA[y,] * WAA)
    SSB[y] <- sum(NAA[y,] * WAA * MatAA) * 0.5
    Catch[y] <- sum(CAA[y,] * WAA)
  }

  # Observation model
  ObsCatch <- Catch * rlnorm(n_yrs, 0, sigma_catch)

  # Survey observations
  SrvIdx <- srv_q * Total_Biom * rlnorm(n_yrs, 0, sigma_srv)

  # Age composition data (simplified multinomial)
  ObsFishAges <- array(0, dim=c(n_yrs, n_ages))
  ObsSrvAges <- array(0, dim=c(n_yrs, n_ages))

  for(y in 1:n_yrs) {
    # Fishery ages
    ObsFishAges[y,] <- rmultinom(1, fish_ISS, CAA[y,])
    # Survey ages
    ObsSrvAges[y,] <- rmultinom(1, srv_ISS, NAA[y,] * srv_sel)
  }

  # Return list of true and observed values
  return(list(
    NAA = NAA,
    CAA = CAA,
    FAA = FAA,
    SSB = SSB,
    Total_Biom = Total_Biom,
    Catch = Catch,
    ObsCatch = ObsCatch,
    SrvIdx = SrvIdx,
    ObsFishAges = ObsFishAges,
    ObsSrvAges = ObsSrvAges,
    fish_sel = fish_sel,
    srv_sel = srv_sel,
    WAA = WAA,
    MatAA = MatAA,
    M = M,
    srv_q = srv_q,
    rec_devs = rec_devs,
    init_devs = init_devs
  ))
}

# Set up simulation -------------------------------------------------------------
years <- 1:20
ages <- 1:15
WAA <- 2 / (1 + exp(-0.8 * (ages - 3)))
MatAA <- 1 / (1 + exp(-1 * (ages - 5)))
sigma_R <- 0.3
sigma_Catch <- 0.001
sigma_SrvIdx <- 0.3
Fmort <- c(seq(0.02, 0.3, length.out = 10), seq(0.3, 0.05, length.out = 10))

# First, simulate some data for the model
sim <- sim_pop_model(years = years,
                     ages = ages,
                     WAA = WAA,
                     MatAA = MatAA,
                     mean_Rec = 50,
                     sigma_R = sigma_R,
                     sigma_catch = sigma_Catch,
                     sigma_srv = sigma_SrvIdx,
                     fish_ISS = 1e5,
                     srv_ISS = 1e5,
                     M = 0.3,
                     fish_sel = 1 / (1 + exp(-2.5 * (ages - 6))),
                     srv_sel = 1 / (1 + exp(-2 * (ages - 3))),
                     Fmort = Fmort,
                     srv_q = 1)


# Set up Rceattle data -------------------------------------------------------------
library(Rceattle)
data("GOAcod")
simData <- GOAcod

# * Data controls ----
simData$nspp <- 1
simData$styr <- 1
simData$endyr <- 20
simData$projyr <- 30
simData$nsex <- 1
simData$nages <- 15
simData$minage <- 1
simData$nlengths <- 15
simData$pop_wt_index <- 1
simData$ssb_wt_index <- 1
simData$pop_age_transition_index <- 1

# * Fleet control ----
simData$fleet_control <- simData$fleet_control[c(1,3),] # BT and Trawl Fishery are both simple logistic
simData$fleet_control$Fleet_name <- c("Survey", "Fishery")
simData$fleet_control$Fleet_code <- 1:2
simData$fleet_control$Selectivity_index <- 1:2
simData$fleet_control$Weight_index <- 1

# * Index data ----
simData$index_data <- data.frame(Fleet_name = "Survey",
                                 Fleet_code = 1,
                                 Species = 1,
                                 Year = 1:20,
                                 Month = 0,
                                 Selectivity_block = 1,
                                 Q_block = 1,
                                 Observation = sim$SrvIdx,
                                 Log_sd = sigma_SrvIdx)

# * Catch data ----
simData$catch_data <- data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Year = 1:20,
                                 Month = 0,
                                 Selectivity_block = 1,
                                 Catch = sim$ObsCatch,
                                 Log_sd = sigma_Catch)

# * Comp data ----
# - Index
tmp <- sim$ObsSrvAges
colnames(tmp) <- paste0("Comp_",1:15)
index_comp <- cbind(data.frame(Fleet_name = "Survey",
                               Fleet_code = 1,
                               Species = 1,
                               Sex = 0,
                               Age0_Length1 = 0,
                               Year = 1:20,
                               Month = 0,
                               Sample_size = rowSums(tmp)),
                    tmp
)

# - Fishery
tmp <- sim$ObsFishAges
colnames(tmp) <- paste0("Comp_",1:15)
fishery_comp <- cbind(data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Sex = 0,
                                 Age0_Length1 = 0,
                                 Year = 1:20,
                                 Month = 0,
                                 Sample_size = rowSums(tmp)),
                      tmp
)

simData$comp_data <- rbind(index_comp, fishery_comp)

# * Empirical selectivity ----
simData$emp_sel[] <- NA

# * Fixed numbers ----
simData$NByageFixed[] <- NA


# * Age transition matrix ----
tmp <- as.data.frame(diag(1,15))
colnames(tmp) <- paste0("Length_",1:15)
simData$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                             Age_transition_index = 1,
                                             Species = 1,
                                             Sex = 0,
                                             Age = 1:15),
                                  tmp
)


# * Age error ----
tmp <- as.data.frame(diag(1,15))
colnames(tmp) <- paste0("Obs_age",1:15)
simData$age_error <- cbind(data.frame(Species = 1,
                                      True_age = 1:15),
                           tmp
)

# * Weight-at-age ----
WAA <- as.data.frame(matrix(WAA, ncol = 15))
colnames(WAA) <- paste0("Age",1:15)
simData$wt <- cbind(data.frame(Wt_name = "Base",
                               Wt_index = 1,
                               Species = 1,
                               Sex = 0,
                               Year = 0),
                    WAA
)

# * Maturity ----
MatAA <- as.data.frame(matrix(MatAA, ncol = 15))
colnames(MatAA) <- paste0("Age",1:15)
simData$pmature <- cbind(data.frame(Species = 1),
                         MatAA
)

# * Sex ratio ----
sexratio <- as.data.frame(matrix(0.5, ncol = 15))
colnames(sexratio) <- paste0("Age",1:15)
simData$sex_ratio <- cbind(data.frame(Species = 1),
                           sexratio
)

# * Mortality ----
mort <- as.data.frame(matrix(0.3, ncol = 15))
colnames(mort) <- paste0("Age",1:15)
simData$M1_base <- cbind(data.frame(Species = 1,
                                    Sex = 0),
                         mort
)

# * Environmental data ----
simData$env_data <- data.frame(Year = 1:20,
                               EnvData = 1)


# * Relative foraging rate (days) ----
pyrs <- as.data.frame(matrix(1, nrow = 20, ncol = 15))
colnames(pyrs) <- paste0("Age",1:15)
simData$Pyrs <- cbind(data.frame(Species = 1,
                                 Sex = 0,
                                 Year = 1:20),
                      pyrs
)


# Fit Rceattle -------------------------------------------------------------
ss_run <- Rceattle::fit_mod(data_list = simData,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

plot(x = sim$SSB, y = ss_run$quantities$ssb[1,1:20]); abline(1,1)
plot(x = sim$Total_Biom, y = ss_run$quantities$biomass[1,1:20], ylab = "Rceattle biomass", xlab = "True biomass"); abline(1,1)
plot(x = sim$NAA[1:20,1], y = ss_run$quantities$R[1,1:20]); abline(1,1)
