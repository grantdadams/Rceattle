library(Rceattle)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
CSL_SBF_ATF_hakedata <- Rceattle::read_data(file = "hake_yr24_SBF_ATF_CSL_Final_projections.xlsx")
CSL_SBF_ATF_hakedata$projyr <- 2066
CSL_SBF_ATF_hakedata$fleet_control$proj_F_prop # only using for hake fishery

# Check current weight data structure
head(CSL_SBF_ATF_hakedata$weight)

# Find which indices use year = 0 (time-invariant)
CSL_SBF_ATF_hakedata$weight %>%
  filter(Year == 0) %>%
  select(Wt_name, Wt_index, Sex, Year) %>%
  distinct()

# Expand year = 0 rows to cover all hindcast years
hindcast_years <- 1966:2024

# Get the time-invariant rows
wt_invariant <- CSL_SBF_ATF_hakedata$weight %>%
  filter(Year == 0)

# Get the time-varying rows (hake)
wt_varying <- CSL_SBF_ATF_hakedata$weight %>%
  filter(Year != 0)

# Expand time-invariant rows across all hindcast years
wt_expanded <- map_dfr(hindcast_years, function(yr) {
  wt_invariant %>% mutate(Year = yr)
})

# Combine back together
CSL_SBF_ATF_hakedata_mse <- CSL_SBF_ATF_hakedata
CSL_SBF_ATF_hakedata_mse$weight <- bind_rows(
  wt_varying,
  wt_expanded
) %>%
  arrange(Wt_index, Sex, Year)

# Verify all indices now span full hindcast period
CSL_SBF_ATF_hakedata_mse$weight %>%
  group_by(Wt_index, Sex) %>%
  summarise(
    min_year = min(Year),
    max_year = max(Year),
    n_years  = n_distinct(Year)
  ) %>%
  print()

# All should show styr to endyr
CSL_SBF_ATF_hakedata<- CSL_SBF_ATF_hakedata_mse


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Operating Models ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Single-species with fixed M
ss_run <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = TRUE,
                            verbose = 1)

ss_run$quantities$jnll # 2469.184
ss_run$quantities$M1 #0.1990000

#  Single-species with estimated M
ss_run_M <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                              inits = ss_run$estimated_params, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 1, # Estimate hindcast only
                              random_rec = FALSE, # No random recruitment
                              M1Fun = build_M1(M1_model = 1),
                              msmMode = 0, # Single species mode
                              phase = TRUE,
                              verbose = 1)

ss_run_M$quantities$jnll # 2465.417
ss_run_M$quantities$M1 #0.2658643

# Multispecies model
ms_run <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                            inits = ss_run_M$estimated_params, # Initial parameters from single species ests
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate hindcast only
                            niter = 3, # 10 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            M1Fun = build_M1(updateM1 = TRUE), # Fix residual M to values in data
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

ms_run$quantities$jnll # 2492.357
ms_run$quantities$M1 #0.1990000

ms_run_M <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                            inits = ss_run_M$estimated_params, # Initial parameters from single species ests
                            M1Fun = build_M1(M1_model = 1,  #estimate mortality!
                                             updateM1 = FALSE,
                                             M1_use_prior = FALSE,
                                             M2_use_prior = FALSE),
                            file = NULL, # Don't save
                            estimateMode = 1, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            verbose = 1)

ms_run_M$quantities$jnll # 2549.51
ms_run_M$quantities$M1 #0.3024391

# Plot OMs:
mod_list <- list(ss_run, ss_run_M, ms_run, ms_run_M)
mod_names <- c("SS", "SS-M", "MS", "MS-M")
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)

#############################################################################################################
## LOGNORMAL SUITABILITY MODEL + DM LIKELIHOOD DISTRIBUTION
#############################################################################################################

# Create initial parameter list:
test_data <- CSL_SBF_ATF_hakedata

inits = ms_run_M$estimated_params
map = ms_run_M$map # gam_a, gam_b, and log_phi are turned off here

# Create a list prey size preference
#ORIGINAL (TRUE) VALUES
inits$log_gam_a = c(0, 3.7, 3.1, 0)  # Mean log weight ratio for ATF, 0 for other species (pred/prey)
inits$log_gam_b = c(0, 1.83, 1.120, 0)

# Set vulnerability matrix
inits$log_phi #Currently all set to 0.5 (keep it)
inits$log_phi[1,2] <- -999 #Fixing so hake do not prey on ATF
inits$log_phi[1,3] <- -999 #Fixing so hake do not prey on SBF
inits$log_phi[1,4] <- -999 #Fixing so hake do not prey on CSL

inits$log_phi[2,2] <- -999 # Set ATF do not feed on ATF
inits$log_phi[2,3] <- -999 # Set ATF do not feed on SBF
inits$log_phi[2,4] <- -999 # Set ATF do not feed on CSL

inits$log_phi[3,3] <- -999 # Set SBF do not feed on SBF
inits$log_phi[3,2] <- -999 # Set SBF do not feed on ATF
inits$log_phi[3,4] <- -999 # Set SBF do not feed on CSL

inits$log_phi[4,4] <- -999 # Set CSL do not feed on CSL
inits$log_phi[4,2] <- -999 # Set CSL do not feed on ATF
inits$log_phi[4,3] <- -999 # Set CSL do not feed on SBF


# Do this to estimate vulnerability and log_phi :
map$mapList$log_phi[] <- 1:length(map$mapList$log_phi) # Unique number for each parameter
map$mapList$log_phi[1,1] <- NA #so we dont estimate on hake on hake
map$mapList$log_phi[1,2] <- NA #so we dont estimate on hake on ATF
map$mapList$log_phi[1,3] <- NA #so we dont estimate on hake on SBF
map$mapList$log_phi[1,4] <- NA #so we dont estimate on hake on CSL

map$mapList$log_phi[2,2] <- NA #so we dont estimate ATF on ATF
map$mapList$log_phi[2,3] <- NA #so we dont estimate on SBF on SBF
map$mapList$log_phi[2,4] <- NA #so we dont estimate on ATF on CSL

map$mapList$log_phi[3,2] <- NA #so we dont estimate SBF on ATF
map$mapList$log_phi[3,3] <- NA #so we dont estimate SBF on SBF
map$mapList$log_phi[3,4] <- NA #so we dont estimate SBF on CSL

map$mapList$log_phi[4,1] <- NA #so we dont estimate CSL on hake
map$mapList$log_phi[4,2] <- NA #so we dont estimate CSL on ATF
map$mapList$log_phi[4,3] <- NA #so we dont estimate CSL on SBF
map$mapList$log_phi[4,4] <- NA #so we dont estimate CSL on CSL

map$mapFactor$log_phi <- factor(map$mapList$log_phi)

# RUN MODEL
run_ms_CSL_Mest_prior<- Rceattle::fit_mod(data_list = test_data,
                                          inits = inits, # Initial parameters from single species ests
                                          map = map,
                                          M1Fun = build_M1(M1_model = 1,
                                                           #updateM1 = TRUE,
                                                           M1_use_prior = TRUE,
                                                           M_prior = 0.2,
                                                           M_prior_sd = 0.1),
                                          file = NULL, # Don't save
                                          estimateMode = 1, # estimate hindcast only
                                          niter = 3, # 3 iterations around population and predation dynamics
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = 1, # MSVPA based
                                          loopnum = 5,
                                          phase = TRUE,
                                          suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                          initMode = 2,
                                          verbose = 1)

run_ms_CSL_Mest_prior$quantities$jnll #2769.784
run_ms_CSL_Mest_prior$quantities$M1 #0.2729866

# Plot OMs:
mod_list <- list(ss_run, ss_run_M, ms_run, ms_run_M, run_ms_CSL_Mest_prior)
mod_names <- c("SS", "SS-M", "MS", "MS-M", "MS-M-LN")
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Multi-species harvest control rules ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# -- F that acheives 40% of SB0, where SB0 is derived from projecting all species simultaneously under no fishing
ms_run_fb40_M <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_prior$data_list,
                                          inits = run_ms_CSL_Mest_prior$estimated_params, # Initial parameters from  multispecies ests
                                          map = run_ms_CSL_Mest_prior$map,
                                          M1Fun = build_M1(M1_model = 1,
                                                           #updateM1 = TRUE,
                                                           M1_use_prior = TRUE,
                                                           M_prior = 0.2,
                                                           M_prior_sd = 0.1),
                                          file = NULL, # Don't save
                                          estimateMode = 0, # estimate
                                          niter = 3, # 3 iterations around population and predation dynamics
                                          random_rec = FALSE, # No random recruitment
                                          msmMode = 1, # MSVPA based
                                          loopnum = 5,
                                          phase = TRUE,
                                          suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                          initMode = 2,
                                          HCR = build_hcr(HCR = 3, # Constant F HCR
                                                 DynamicHCR = FALSE, # Use dynamic reference points
                                                 Ftarget = 0.4), # F that achieves 40% SB0
                                          verbose = 1)

ms_run_fb40_M$quantities$jnll #2769.784
ms_run_fb40_M$quantities$M1 #0.2729866

# 2. Dynamic BRPs
ms_run_fb40_M_dynamic <- Rceattle::fit_mod(data_list = run_ms_CSL_Mest_prior$data_list,
                                   inits = run_ms_CSL_Mest_prior$estimated_params, # Initial parameters from  multispecies ests
                                   map = run_ms_CSL_Mest_prior$map,
                                   M1Fun = build_M1(M1_model = 1,
                                                    #updateM1 = TRUE,
                                                    M1_use_prior = TRUE,
                                                    M_prior = 0.2,
                                                    M_prior_sd = 0.1),
                                   file = NULL, # Don't save
                                   estimateMode = 0, # estimate
                                   niter = 3, # 3 iterations around population and predation dynamics
                                   random_rec = FALSE, # No random recruitment
                                   msmMode = 1, # MSVPA based
                                   loopnum = 5,
                                   phase = TRUE,
                                   suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                   initMode = 2,
                                   HCR = build_hcr(HCR = 3, # Constant F HCR
                                                   DynamicHCR = TRUE, # Use dynamic reference points
                                                   Ftarget = 0.4), # F that achieves 40% SB0
                                   verbose = 1)

ms_run_fb40_M_dynamic$quantities$jnll #2769.784
ms_run_fb40_M_dynamic$quantities$M1 #0.2729866

# -- Plot
mod_list <- list(ms_run_fb40_M, ms_run_fb40_M_dynamic)
model_names <- c("F40","F40dyn")
plot_biomass(mod_list, model_names = model_names, incl_proj = TRUE)
plot_ssb(mod_list, model_names = model_names, incl_proj = TRUE)
plot_recruitment(mod_list, model_names = model_names, incl_proj = TRUE)
plot_catch(mod_list, model_names = model_names, incl_proj = TRUE)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Single species harvest control rules ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# -- Constant F as a percentage of SB0
ss_run_fb0 <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                                inits = ss_run$estimated_params, # Initial parameters from ss_run
                                estimateMode = 2, # Run projection only
                                HCR = build_hcr(HCR = 3, # Constant F HCR
                                                DynamicHCR = FALSE, # Use dynamic reference points
                                                Ftarget = 0.4), # F that achieves 40% SB0
                                msmMode = 0, # Single species mode
                                phase = TRUE,
                                verbose = 1)

ss_run_dynamicfb0 <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                                       inits = ss_run$estimated_params, # Initial parameters from ss_run
                                       estimateMode = 2, # Run projection only
                                       HCR = build_hcr(HCR = 3, # Constant F HCR
                                                       DynamicHCR = TRUE, # Use dynamic reference points
                                                       Ftarget = 0.4), # F that achieves 40% SB0
                                       msmMode = 0, # Single species mode
                                       verbose = 1)


# -- PFMC Category 1
ss_run_Cat1 <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                                 inits = ss_run$estimated_params, # Initial parameters from ss_run
                                 estimateMode = 2, # Run projection only
                                 HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                 Flimit = 0.45, # F45%
                                                 Ptarget = 0.4, # Target is 40% B0
                                                 Plimit = 0.1, # No fishing when SB<SB10
                                                 Pstar = 0.45,
                                                 Sigma = 0.5),
                                 msmMode = 0, # Single species mode
                                 verbose = 1)

ss_run_dynamicCat1 <- Rceattle::fit_mod(data_list = CSL_SBF_ATF_hakedata,
                                        inits = ss_run$estimated_params, # Initial parameters from ss_run
                                        estimateMode = 2, # Run projection only
                                        HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                        DynamicHCR = TRUE, # Use dynamic reference points
                                                        Flimit = 0.45, # F45%
                                                        Ptarget = 0.4, # Target is 40% SB0
                                                        Plimit = 0.1, # No fishing when SB<SB10
                                                        Pstar = 0.45,
                                                        Sigma = 0.5),
                                        msmMode = 0, # Single species mode
                                        verbose = 1)


# -- Plot
mod_list <- list(ss_run, ss_run_fb0, ss_run_Cat1)
model_names <- c("F=0","F 40% B0", "PFMC Cat 1")
plot_biomass(mod_list, model_names = model_names, incl_proj = TRUE, species =1, add_ci = TRUE)
plot_ssb(mod_list, model_names = model_names, incl_proj = TRUE, species =1, add_ci = TRUE)
plot_depletionSSB(mod_list, model_names = model_names, incl_proj = TRUE, species =1)

dynamic_mod_list <- list(ss_run, ss_run_dynamicfb0, ss_run_dynamicCat1)
dynamic_model_names <- c("F=0","F 40% B0", "Fspr 40%", "PFMC Cat 1")
plot_biomass(dynamic_mod_list, model_names = dynamic_model_names, incl_proj = TRUE, species = 1)
plot_ssb(dynamic_mod_list, model_names = dynamic_model_names, incl_proj = TRUE, species = 1)
plot_depletionSSB(mod_list, model_names = model_names, incl_proj = TRUE, species =1)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Management strategy evaluation ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# -- No F
# - MS-OM: SS-EM No F
mse1 <- run_mse(om = run_ms_CSL_Mest_prior, em = ss_run, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, sample_rec = TRUE, cap = c(1500000))

# - SS-OM: SS-EM No F
mse2 <- run_mse(om = ss_run_M, em = ss_run, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, sample_rec = TRUE, cap = c(1500000))

##Check plots
plot_ssb(mse2$Sim_1$OM, species =1)
plot_ssb(mse2$Sim_1$OM_no_F, species =1)
plot_ssb(mse2$Sim_1$EM$EM, species =1)

mse_summary(mse1)
mse_summary(mse2)

# -- PFMC Category 1 HCRs
# - MS-OM: SS-EM Tier 3 HCR
mse3 <- run_mse(om = ms_run, em = ss_run_Cat1, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, sample_rec = TRUE, cap = c(1500000))

# - SS-OM: SS-EM Tier 3 HCR
mse4 <- run_mse(om = ss_run_M, em = ss_run_Cat1, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, sample_rec = TRUE, cap = c(1500000))

# - MS-OM: SS-EM dynamic Tier 3 HCR
mse5 <- run_mse(om = ms_run, em = ss_run_dynamicCat1, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, sample_rec = TRUE, cap = c(1500000))

# - SS-OM: SS-EM dynamic Tier 3 HCR
mse6 <- run_mse(om = ss_run_M, em = ss_run_dynamicCat1, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, sample_rec = TRUE, cap = c(1500000))






# 1. Average catch
avg_catch <- mean(catch_proj, na.rm = TRUE)

# 2. Catch interannual variability (IAV)
iav <- sqrt(mean(diff(catch_proj)^2)) / mean(catch_proj)

# 3. Terminal depletion relative to B0
terminal_depletion <- ssb_proj[length(ssb_proj)] /
  B0_proj[length(B0_proj)]

# 4. Probability SSB < 10% B0 (overfished)
p_overfished <- mean(ssb_proj < 0.1 * B0_proj, na.rm = TRUE)

# 5. Probability SSB < 40% B0 (below target)
p_below_target <- mean(ssb_proj < 0.4 * B0_proj, na.rm = TRUE)

# Summary
cat("=== MSE Performance Metrics ===\n")
cat("Average catch:        ", round(avg_catch / 1e6, 3), "million t\n")
cat("Catch IAV:            ", round(iav, 3), "\n")
cat("Terminal depletion:   ", round(terminal_depletion, 3), "\n")
cat("P(overfished):        ", round(p_overfished, 3), "\n")
cat("P(below target):      ", round(p_below_target, 3), "\n")


###### TETSING ####
# Check SSB trajectory during projection
SSB_proj <- hake_ms_F0$quantities$ssb

# Dims: [species, sex, year]
dim(SSB_proj)

# Build dataframe for hake (species 1, sex 1)
years_full <- run_ms_CSL_Mest_prior$data_list$styr:hake_data_F0$projyr
n_years    <- length(years_full)

SSB_df <- data.frame(
  Year = years_full,
  SSB  = SSB_proj[1, 1:n_years]
) %>%
  mutate(
    SSB_lag    = lag(SSB),
    pct_change = abs(SSB - SSB_lag) / SSB_lag * 100,
    Era        = ifelse(Year <= run_ms_CSL_Mest_test_W$data_list$endyr,
                        "Hindcast", "Projection")
  )

# Plot to see if equilibrium is reached within 50 years
ggplot(SSB_df, aes(x = Year, y = SSB / 1e6, color = Era)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 2019, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Hindcast" = "#0072B2",
                                "Projection" = "#D55E00")) +
  theme_minimal(base_size = 14) +
  labs(
    y     = "Spawning biomass (million t)",
    x     = "Year",
    color = NULL,
    title = "Hake SSB under F=0 — finding multispecies B0"
  ) +
  theme(legend.position = "bottom")

# Find equilibrium year and B0
equil_year <- SSB_df %>%
  filter(Era == "Projection", pct_change < 0.1) %>%
  slice(1) %>%
  pull(Year)

B0_ms <- SSB_df %>%
  filter(Year == equil_year) %>%
  pull(SSB)

cat("Equilibrium reached:", equil_year, "\n")
cat("Multispecies B0:    ", round(B0_ms / 1e6, 3), "million t\n")




#########################################
# EM has two levels:
# mse2$Sim_1$EM$EM          = the initial EM fit to hindcast data
# mse2$Sim_1$EM$`OM_Sim_1. EM_yr_XXXX` = annual EM updates during projection

# ── Initial EM (hindcast fit) ──────────────────────────
names(mse2$Sim_1$EM$EM$quantities)  # should have full quantities list

# ── Annual EM updates (projection period) ─────────────
# Each year has a shorter quantities list (24 items vs 77)
names(mse2$Sim_1$EM$`OM_Sim_1. EM_yr_2024`$quantities)

# ── Extract SSB from OM and EM for comparison ──────────
years_full <- CSL_SBF_ATF_hakedata$styr:CSL_SBF_ATF_hakedata$projyr
years_proj <- (CSL_SBF_ATF_hakedata$endyr + 1):CSL_SBF_ATF_hakedata$projyr

# OM true SSB (full time series)
ssb_OM <- data.frame(
  Year  = years_full,
  SSB   = mse2$Sim_1$OM$quantities$ssb[1, ],
  Model = "OM (true)"
)

# OM no F SSB (unfished - for true B0)
ssb_OM_noF <- data.frame(
  Year  = years_full,
  SSB   = mse2$Sim_1$OM_no_F$quantities$ssb[1, ],
  Model = "OM no F (B0)"
)

# EM hindcast SSB
ssb_EM_hindcast <- data.frame(
  Year  = CSL_SBF_ATF_hakedata$styr:CSL_SBF_ATF_hakedata$endyr,
  SSB   = mse2$Sim_1$EM$EM$quantities$ssb[1, ],
  Model = "EM (perceived)"
)

# EM annual SSB estimates during projection
# Each annual EM gives SSB estimate for its terminal year
em_years <- years_proj
em_names <- paste0("OM_Sim_1. EM_yr_", em_years)

ssb_EM_proj <- map2_dfr(em_years, em_names, function(yr, nm) {

  em_yr <- mse2$Sim_1$EM[[nm]]

  if(is.null(em_yr$quantities$ssb)) return(NULL)

  # Get terminal year SSB from each annual assessment
  ssb_vals <- em_yr$quantities$ssb[1, ]
  n <- length(ssb_vals)

  data.frame(
    Year  = yr,
    SSB   = ssb_vals[n],   # terminal year estimate
    Model = "EM (perceived)"
  )
})

# Combine EM SSB across hindcast and projection
ssb_EM_full <- bind_rows(ssb_EM_hindcast, ssb_EM_proj)

# ── Plot OM vs EM SSB ──────────────────────────────────
bind_rows(ssb_OM, ssb_EM_full) %>%
  ggplot(aes(x = Year, y = SSB / 1e6, color = Model)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = CSL_SBF_ATF_hakedata$endyr,
             linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(
    "OM (true)"      = "#D55E00",
    "EM (perceived)" = "#0072B2"
  )) +
  theme_minimal(base_size = 14) +
  labs(y = "SSB (million t)", color = NULL,
       title = "MSE: true vs perceived SSB") +
  theme(legend.position = "bottom")

# ── Extract B0 and depletion ───────────────────────────
# True B0 from OM_no_F
B0_true <- mse2$Sim_1$OM_no_F$quantities$ssb[1, ]

# EM perceived B0 from each annual assessment
B0_EM_proj <- map2_dfr(em_years, em_names, function(yr, nm) {
  em_yr <- mse2$Sim_1$EM[[nm]]
  if(is.null(em_yr$quantities$SB0)) return(NULL)
  data.frame(
    Year = yr,
    B0   = em_yr$quantities$SB0[length(em_yr$quantities$SB0)]
  )
})

# ── Performance metrics ────────────────────────────────
proj_idx  <- which(years_full %in% years_proj)
ssb_proj  <- mse2$Sim_1$OM$quantities$ssb[1, proj_idx]
catch_proj <- mse2$Sim_1$OM$quantities$catch_hat[1, 1, proj_idx]
B0_proj   <- B0_true[proj_idx]

cat("=== Performance metrics ===\n")
cat("Avg catch (million t):  ",
    round(mean(catch_proj, na.rm = TRUE) / 1e6, 3), "\n")
cat("Catch IAV:              ",
    round(sqrt(mean(diff(catch_proj)^2)) /
            mean(catch_proj, na.rm = TRUE), 3), "\n")
cat("Terminal depletion:     ",
    round(ssb_proj[length(ssb_proj)] /
            B0_proj[length(B0_proj)], 3), "\n")
cat("P(SSB < 40% B0):        ",
    round(mean(ssb_proj < 0.4 * B0_proj, na.rm = TRUE), 3), "\n")
cat("P(SSB < 10% B0):        ",
    round(mean(ssb_proj < 0.1 * B0_proj, na.rm = TRUE), 3), "\n")

######################################################################

# Load your fitted multispecies model
# (equivalent to ms_run in Grant's code but your hake model)
# You already have: run_ms_CSL_Mest_test_W

# Set projection end year
hake_data_proj <- run_ms_CSL_Mest_prior$data_list
hake_data_proj$projyr <- 2066
run_ms_CSL_Mest_prior$data_list$endyr

# Check projection F settings for predator fleets
# You want predator fleets to use their SS-projected abundance
# not be estimated during projection
hake_data_proj$fleet_control$proj_F_prop #we are projecting the hake fishery fleet

# Primary run: F=0 on hake, predators follow SS projections
hake_ms_F0 <- Rceattle::fit_mod(data_list = hake_data_proj,
                                inits = run_ms_CSL_Mest_prior$estimated_params, # Initial parameters from single species ests
                                map = run_ms_CSL_Mest_prior$map,
                                M1Fun = build_M1(M1_model = 1,
                                                 #updateM1 = TRUE,
                                                 M1_use_prior = TRUE,
                                                 M_prior = 0.2,
                                                 M_prior_sd = 0.1),
                                file = NULL, # Don't save
                                estimateMode = 0, # only projections
                                niter = 3, # 3 iterations around population and predation dynamics
                                random_rec = FALSE, # No random recruitment
                                msmMode = 1, # MSVPA based
                                loopnum = 5,
                                phase = TRUE,
                                suitMode = c(0, 4, 4, 0), # empirical + LN suitability
                                initMode = 2,
                                HCR = build_hcr(
                                  HCR = 3,               # constant F
                                  DynamicHCR = FALSE,    # static B0 first
                                  Ftarget = 0.4),         # F that achieves 40% SB0
                                verbose = 1)


hake_ms_F0$quantities$jnll

plot_biomass(ss_run, incl_proj = TRUE, species =1)
plot_biomass(hake_ms_F0, incl_proj = TRUE, species =1)

plot_ssb(ss_run, incl_proj = TRUE, species =1)
plot_ssb(hake_ms_F0, incl_proj = TRUE, species =1)

# Find the correct terminal index
terminal_idx <- length(hake_ms_F0$quantities$B0[1,])
cat("Terminal index:", terminal_idx, "\n")
cat("Corresponds to year:",
    run_ms_CSL_Mest_prior$data_list$styr + terminal_idx - 1, "\n")

brp_comparison <- function(model, model_name) {
  # Note: Use index [1] for the first species (Hake)
  df <- data.frame(Value = c(
    (model$quantities$B0[terminal_idx] / 1000000),      # Unfished Biomass (million tons)
    (model$quantities$SB0[terminal_idx] / 1000000),     # Unfished Spawning Biomass
    (model$quantities$R0[1] / 1000000),      # Equilibrium Recruitment
    model$quantities$Ftarget[1],             # HCR prescribed F
    model$quantities$proj_F[1, terminal_idx],           # Projected F for species 1, fleet 1
    model$quantities$SPRlimit[1],            # SPR at Flimit
    model$quantities$SPRtarget[1]            # SPR at Ftarget
  ))

  row.names(df) <- c("B0 (Mt)", "SB0 (Mt)", "R0 (millions)", "Ftarget", "Fproj", "SPRlimit", "SPRtarget")
  colnames(df) <- model_name
  df <- round(df, 3) # Increased to 3 decimal places for F rates
  return(df)
}

# Run the comparison
brp_comparison(model = hake_ms_F0, model_name = "MS")
brp_comparison(model = ss_run, model_name = "MS")


# Extract SSB trajectory during projection
SSB_proj <- hake_ms_F0$quantities$ssb

# Find equilibrium year - where SSB stabilizes
SSB_proj_df <- data.frame(
  Year = hake_data_proj$styr:hake_data_proj$projyr,
  SSB  = SSB_proj[1, 1:87]  # species 1
) %>%
  mutate(
    SSB_lag  = lag(SSB),
    pct_change = abs(SSB - SSB_lag) / SSB_lag * 100
  )

# Plot to visually find equilibrium
ggplot(SSB_proj_df, aes(x = Year, y = SSB)) +
  geom_line(linewidth = 1.2, color = "#D55E00") +
  geom_vline(
    data = SSB_proj_df %>% filter(pct_change < 0.05),
    aes(xintercept = min(Year)),
    linetype = "dashed", color = "grey50"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    y = "Spawning biomass (t)",
    x = "Year",
    title = "Hake SSB under F=0 — finding multispecies B0",
    caption = "Dashed line = equilibrium (annual change < 0.1%)"
  )

# Extract B0 at equilibrium
B0_multispecies <- SSB_proj_df %>%
  filter(pct_change < 0.05) %>%  # converged
  slice(1) %>%
  pull(SSB)

# Compare to single-species B0
dim(ss_run$quantities$B0) #the single spcies model run up to 2100, we need to compare to 2066 value?
B0_ss <- ss_run$quantities$B0[1]
# adjust to your ss model object name

cat("Single-species B0:  ", round(B0_ss / 1e6, 3), "million t\n")
cat("Multispecies B0:    ", round(B0_multispecies / 1e6, 3), "million t\n")
cat("Ratio (MS/SS):      ", round(B0_multispecies / B0_ss, 3), "\n")
cat("Reduction from pred:", round((1 - B0_multispecies/B0_ss)*100, 1), "%\n")



# Check what projyr values exist
run_ms_CSL_Mest_prior$data_list$projyr   # original fitted model
hake_ms_F0$data_list$projyr               # what the F=0 run actually used

# Also check the quantities array dimensions
# This tells you how many years the model actually ran
dim(hake_ms_F0$quantities$ssb)
# 3rd dimension = total years (hindcast + projection)

# Calculate what years that covers
styr   <- hake_ms_F0$data_list$styr
n_yrs  <- dim(hake_ms_F0$quantities$biomassSSB)[3]
endyr_actual <- styr + n_yrs - 1
cat("Model actually ran to:", endyr_actual, "\n")




