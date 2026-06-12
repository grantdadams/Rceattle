# Comparison of Rceattle and WHAM fit to conditional age-at-length (CAAL) data
# Modified from https://giancarlomcorrea.netlify.app/labs/OFI_WK_2023/examples/case3.R
# Note that the plus group in year 1 is wrong on WHAM line 1720 of the cpp and needs adjustment

# OSA residuals: Rceattle vs WHAM --------------------------------------------
# Cross-check Rceattle's one-step-ahead (OSA) residuals against WHAM's (the
# reference implementation that Rceattle's composition decomposition was ported
# from). Both packages are fit to the same GOAcod + whamGrowthData; with the
# fits matched (above), the aggregate index/catch OSA residuals (lognormal) and
# the composition residuals should fall on the 1:1 line. Conditioning order can
# differ slightly between packages, so expect a tight correlation rather than an
# exact equality.
#
# NOTE: manual comparison using the WHAM growth fork installed at the top of this
# script (GiancarloMCorrea/wham, ref = 'growth'). Its OSA residuals come from
# make_osa_residuals(), and its composition decomposition (src/age_comp_osa.hpp)
# is the same conditional binomial / beta-binomial that Rceattle ported into
# comp_osa.hpp -- so the residuals should agree (up to conditioning order).

# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))
library(readxl)
library(wham)
library(ggplot2)
library(dplyr)
source("tests/comparison/WHAM_growth_comparison/helper.R")
dir.create("tests/comparison/WHAM_growth_comparison")
runmodels = FALSE

# WHAM ------------------------------------------------------------------
catch_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 1)
index_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 2)
fsh_lcomp_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 3)
index_lcomp_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 4)
caal_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 5)
maturity_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 6)


# * Load model ----
wham_model <- readRDS(file = "tests/comparison/WHAM_growth_comparison/wham_model.RDS")


# WHAM OSA residuals on the existing fit. Growth-fork `$osa` type labels:
# "logcatch"/"logindex" (aggregate), "catchpal"/"indexpal" (length comp),
# "catchcaal"/"indexcaal" (conditional age-at-length).
wham_osa_fit <- wham::make_osa_residuals(
  wham_model, osa.opts = list(method = "oneStepGaussianOffMode", parallel = TRUE))
wham_osa <- wham_osa_fit$osa            # = input$data$obs + a `residual` column
print(unique(wham_osa$type))            # confirm the labels for your WHAM version


# Rceattle -------------------------------------------------------------
library(Rceattle)
data("GOAcod")
simData <- GOAcod


# * Make input data ----
# * Data controls
simData$spnames <- "WHAM data"
simData$nspp <- 1
years = 1976:2020 # years
nyrs <- length(years)
simData$styr <- 1976
simData$endyr <- 2020
simData$projyr <- 2030
simData$nsex <- 1
simData$nages <- 10
simData$minage <- 1
lengths = seq(from = 2, to = 130, by = 2) # lengths
simData$nlengths <- length(lengths)
simData$pop_wt_index <- 1
simData$ssb_wt_index <- 1
simData$spawn_month <- 0
simData$alpha_wt_len <- 5.56e-06
simData$beta_wt_len <- 3.2
simData$pop_age_transition_index <- 1
simData$sigma_rec_prior <- 0.6
nages <- simData$nages
nlengths <- simData$nlengths

# * Fleet control
simData$fleet_control <- simData$fleet_control[c(1,3),] # BT and Trawl Fishery are both simple logistic
simData$fleet_control$Fleet_name <- c("Survey", "Fishery")
simData$fleet_control$Fleet_code <- 1:2
simData$fleet_control$Selectivity <- "Logistic"                  # Length-based logistic
simData$fleet_control$Selectivity_dimension <- "Length"          # Length-based logistic
simData$fleet_control$Selectivity_index <- 1:2
simData$fleet_control$Weight_index <- 1
simData$fleet_control$Weight1_Numbers2 <- 1
simData$fleet_control$Month <- c(0, 6)
simData$fleet_control$Bin_first_selected <- 1
simData$fleet_control$Catchability <- c(0, NA)
simData$fleet_control$Q_prior <- c(1, NA)

# * Index data
simData$index_data <- data.frame(Fleet_name = "Survey",
                                 Fleet_code = 1,
                                 Species = 1,
                                 Year = years,
                                 Month = index_df$fr_yr,
                                 Selectivity_block = 1,
                                 Q_block = 1,
                                 Observation = index_df$index,
                                 Log_sd = index_df$cv)

# * Catch data
simData$catch_data <- data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Year = years,
                                 Month = 0,
                                 Selectivity_block = 1,
                                 Catch = catch_df$catch,
                                 Log_sd = catch_df$cv)

# * Comp data
# - Index length comp
tmp <- as.matrix(index_lcomp_df[,3:67])
colnames(tmp) <- paste0("Comp_",1:simData$nlengths)
index_comp <- cbind(data.frame(Fleet_name = "Survey",
                               Fleet_code = 1,
                               Species = 1,
                               Sex = 0,
                               Age0_Length1 = 1,
                               Year = years, # Negative turns data off
                               Month = 0,
                               Sample_size = index_lcomp_df$Nsamp),
                    tmp
)

# - Fishery length comp
tmp <- as.matrix(fsh_lcomp_df[,3:67])
colnames(tmp) <- paste0("Comp_",1:simData$nlengths)
fishery_comp <- cbind(data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Sex = 0,
                                 Age0_Length1 = 1,
                                 Year = years, # Negative turns data off
                                 Month = 0,
                                 Sample_size = fsh_lcomp_df$Nsamp),
                      tmp
)

simData$comp_data <- rbind(index_comp, fishery_comp)


# * CAAL data
# - Index
index_caal <- caal_df
colnames(index_caal) <- c("Year", "Length", "Sample_size", paste0("CAAL_",1:simData$nages))
simData$caal_data <- index_caal |>
  dplyr::mutate(Fleet_name = "Survey",
                Fleet_code = 1,
                Species = 1,
                Sex = 0) |>
  select(Fleet_name, Fleet_code, Species, Sex, Year, Length, Sample_size, paste0("CAAL_",1:simData$nages)) |>
  as.data.frame()
# - No fishery CAAL

# * Empirical selectivity
simData$emp_sel[] <- NA

# * Fixed numbers
simData$NByageFixed[] <- NA


# * Age transition matrix
tmp <- as.data.frame(diag(1,simData$nlengths))[1:simData$nages,]
colnames(tmp) <- paste0("Length_",1:simData$nlengths)
simData$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                             Age_transition_index = 1,
                                             Species = 1,
                                             Sex = 0,
                                             Age = 1:simData$nages),
                                  tmp
)


# * Age error
tmp <- as.data.frame(diag(1,simData$nages))
colnames(tmp) <- paste0("Obs_age",1:simData$nages)
simData$age_error <- cbind(data.frame(Species = 1,
                                      True_age = 1:simData$nages),
                           tmp
)

# * Empirical weight-at-age
WAA <- as.data.frame(matrix(1:simData$nages, ncol = simData$nages))
colnames(WAA) <- paste0("Age",1:simData$nages)
simData$weight <- cbind(data.frame(Wt_name = "Base",
                                   Wt_index = 1,
                                   Species = 1,
                                   Sex = 0,
                                   Year = 0), # 0 fills all years
                        WAA
)

# * Maturity (time-invariant in Rceattle)
MatAA <- matrix(maturity_df[1,2:11], ncol = simData$nages)
colnames(MatAA) <- paste0("Age",1:simData$nages)
simData$maturity <- cbind(data.frame(Species = 1),
                          MatAA
) |>
  as.data.frame()

# * Sex ratio
sexratio <- as.data.frame(matrix(0.5, ncol = simData$nages))
colnames(sexratio) <- paste0("Age",1:simData$nages)
simData$sex_ratio <- cbind(data.frame(Species = 1),
                           sexratio
)

# * Input mortality-at-age
mort <- as.data.frame(matrix(0.35, ncol = simData$nages))
colnames(mort) <- paste0("Age",1:simData$nages)
simData$M1_base <- cbind(data.frame(Species = 1,
                                    Sex = 0),
                         mort
)

# * Environmental data
simData$env_data <- data.frame(Year = years,
                               EnvData = 1)


# * Ration data
pyrs <- as.data.frame(matrix(1, nrow = 1, ncol = simData$nages))
colnames(pyrs) <- paste0("Age",1:simData$nages)
simData$ration_data <- cbind(data.frame(Species = 1,
                                        Sex = 0,
                                        Year = 0), # Zero fills out all years
                             pyrs
)


# * Null Rceattle ----
ss_fix <- Rceattle::fit_mod(data_list = simData,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 3, # Don't estimate
                            growthFun = build_growth(fun = "vonBertalanffy"), # Von Bert
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = FALSE,
                            verbose = 1)


# * Fix-parameters and check derived quantities ----
inits <- ss_fix$estimated_params # get paramter object from NULL model
names(wham_model)
names(wham_model$parList)

# - Rec
inits$rec_pars[1,1] <- wham_model$parList$mean_rec_pars
inits$rec_dev[1,1] <- wham_model$parList$log_N1_pars[1] -  wham_model$parList$mean_rec_pars
inits$rec_dev[1,2:nyrs] <- wham_model$parList$log_NAA[,1] -  wham_model$parList$mean_rec_pars
inits$init_dev[1,] <- wham_model$parList$log_N1_pars[1] -  wham_model$parList$mean_rec_pars # WHAM assumes rec-dev in year 1 is applied to year-1 to year - nages

# - F (random walk in WHAM)
inits$log_F[2,1]  <- wham_model$parList$log_F1
for(y in 2:nyrs){
  inits$log_F[2,y] <- inits$log_F[2,y-1] + wham_model$parList$F_devs[y-1,1]
}
inits$log_Finit[1] <- -Inf

# - Selectivity
selpars <-  wham_model$env$data$selpars_lower + ( wham_model$env$data$selpars_upper -  wham_model$env$data$selpars_lower) / (1.0 + exp(-(wham_model$parList$logit_selpars)))
inits$log_sel_slp[1,,1] <- rev(log(1/selpars[,24]))
inits$sel_inf[1,,1] <- rev(selpars[,23])

# - Q
inits$index_log_q[1] <- log(wham_model$env$data$q_lower + (wham_model$env$data$q_upper - wham_model$env$data$q_lower) / (1 + exp(-wham_model$parList$logit_q)))

# - Growth
inits$log_growth_pars[1,1,1:3] <- wham_model$parList$growth_a[c(1,3,2),1]
inits$growth_log_sd[1,1,] <- wham_model$parList$SDgrowth_par



# Rceattle OSA residuals: estimateMode = 1 (objective differentiable, hindcast
# only) and fit_control(osa = TRUE) so the composition OSA data is built.
ss_osa <- Rceattle::fit_mod(data_list = simData, inits = inits,
                            estimateMode = 1,
                            growthFun = build_growth(fun = "vonBertalanffy"),
                            random_rec = FALSE, msmMode = 0,
                            initMode = "NonEquilibrium",
                            fit_control = fit_control(phase = TRUE, verbose = 0,
                                                      osa = TRUE))
rce_osa <- Rceattle::osa_residuals(ss_osa, source = c("index", "catch", "comp", "caal"))


# Compare ----
# * Aggregate index & catch (one lognormal residual per year) ----
wham_agg_type <- c(catch = "logcatch", index = "logindex")
par(mfrow = c(1, 2))
for (src in c("index", "catch")) {
  rce <- rce_osa[rce_osa$type == src, c("year", "residual")]
  wh  <- wham_osa[wham_osa$type == wham_agg_type[[src]], c("year", "residual")]
  wh$year <- wh$year + simData$styr - 1
  m   <- merge(rce, wh, by = "year", suffixes = c(".rce", ".wham"))
  if (nrow(m) == 0) next
  plot(m$residual.wham, m$residual.rce,
       xlab = paste("WHAM", src, "OSA residual"),
       ylab = paste("Rceattle", src, "OSA residual"),
       main = sprintf("%s OSA  (r = %.3f)", src,
                      stats::cor(m$residual.wham, m$residual.rce)))
  abline(0, 1)
}

# * Composition: ----
# this is a length/growth model, so Rceattle "comp" is length
# composition (WHAM "indexpal"/"catchpal") and "caal" is conditional
# age-at-length (WHAM "indexcaal"/"catchcaal"). WHAM stores the bin in its `age`
# column; match by (year, bin) -- check colnames(wham_osa) for your version:
rce <- rce_osa[rce_osa$type == "comp", c("year", "age_or_length_bin", "residual")]
wh  <- wham_osa[wham_osa$type == "indexpal", c("year", "bin", "residual") ]
wh$year <- wh$year + simData$styr - 1

# Convert Rceattle's bin index to the length value WHAM uses
rce <- as.data.frame(rce)
rce$bin <- 2 * rce$age_or_length_bin
wh$bin  <- as.numeric(wh$bin)

m <- merge(rce, wh, by = c("year", "bin"))
plot(m$residual.y, m$residual.x); abline(0, 1)   # WHAM (x) vs Rceattle (y)

