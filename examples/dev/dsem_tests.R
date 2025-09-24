

library(Rceattle)
library(dplyr)
library(readxl)

# Pollock
data("GOApollock")
# Differences between CEATTLE and SAFE: the SAFE model penalizes the first 7 and last recruitment deviates, while I penalize them all. Besides that its parameterized the exact same.
GOApollock$styr = 1970
m0 <- Rceattle::fit_mod(
  data_list = GOApollock,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # 0 = Estimate, 1 = Dont run estimation
  random_rec = TRUE, # No random recruitment
  msmMode = 0, # 0 = Single species mode, 1 = MSVPA multi-species mode
  verbose = 1, # Silence optimization output
  phase = FALSE) # Use default phasing

years <- dimnames(m0$initial_params$rec_dev)[[2]]
c(est=as.list(m0$sdrep, "Est", report=TRUE)$sigmaR, sd=as.list(m0$sdrep, "Std", report=TRUE)$sigmaR)
c(est=as.list(m0$sdrep, "Est")$R_ln_sd, sd=as.list(m0$sdrep, "Std")$R_ln_sd)
# 1.3715813 0.1418742

ssb0 <- m0$quantities$ssb
sum(ssb0) # 28040230
rec0 <- m0$quantities$R_hat
sum(rec0) # 88547344
plot(years, ssb0[1,])

## Load ESP data
ESPdata <- read.csv('examples/ESP_2024_4SCEAM.csv')
## set number of years for projection

ny_proj <- as.numeric(max(years))-2023

## get them ready
ESPdata_MS  <- ESPdata %>%
  # adding NA before (match assessment begining)
  bind_rows(data.frame(Year=c(1970:1976))) %>% arrange(Year) %>%
  #  sdding NA after (forecast)
  add_row(Sum_OffYOY_Shelikof=rep(NA,ny_proj)) %>%
  # selecting time series of interest
  select(Spr_SST,Wind_NS,Sum_LCopepod_Shelikof,Sum_Euph_Kodiak,
         Spr_Larvae_Shelikof,Sum_OffYOY_Shelikof,Sum_NearYOY_Kodiak,
         Sum_OffYOY_Cond_Shelikof,
         Sum_Juv_EuphDiet,Fal_Adult_Cond_Fishery) %>%
  #scale and log data
  mutate(across(.cols=c(Spr_Larvae_Shelikof,Sum_OffYOY_Shelikof,Sum_NearYOY_Kodiak),.fns=log)) %>%
  mutate(across(.fns =scale,.cols=-c(Spr_Larvae_Shelikof,Sum_OffYOY_Shelikof,Sum_NearYOY_Kodiak)))
# final data set for modeling
ESPdata_MSr <- ESPdata_MS %>% mutate(recdevs=NA) %>% relocate(recdevs)


sem_MapMod = "
  # link, lag, param_name, start_value
  #internal data relation
  Spr_SST -> Spr_SST,1,AR_SST
  Wind_NS -> Wind_NS,1,AR_Wind
  Sum_LCopepod_Shelikof -> Sum_LCopepod_Shelikof,1, AR_Cope
  Sum_Euph_Kodiak -> Sum_Euph_Kodiak,1,AR_Euph
  Spr_Larvae_Shelikof - > Spr_Larvae_Shelikof,1,AR_larvae
  Sum_OffYOY_Cond_Shelikof -> Sum_OffYOY_Cond_Shelikof,1,AR_CondOffYOY
  Sum_OffYOY_Shelikof -> Sum_OffYOY_Shelikof,1,AR_OffYOY
  Sum_NearYOY_Kodiak -> Sum_NearYOY_Kodiak,1,AR_NearYOY
  Sum_Juv_EuphDiet -> Sum_Juv_EuphDiet,1,AR_JuvEuphDiet
  Fal_Adult_Cond_Fishery -> Fal_Adult_Cond_Fishery,1,AR_CondAd
  recdevs <-> recdevs, 0, sigmaR, 1

  #causal relation
  Sum_Euph_Kodiak -> Fal_Adult_Cond_Fishery, 0, Euph_to_CondAd
  Fal_Adult_Cond_Fishery -> Spr_Larvae_Shelikof,1,CondAd_to_Larvae
  #Spr_SST -> Spr_Larvae_Shelikof,0,SST_to_Larvae
  Wind_NS -> Spr_Larvae_Shelikof,0,Wind_to_Larvae
  Spr_Larvae_Shelikof -> Sum_OffYOY_Shelikof,0,Larvae_to_OffYOY
  Sum_OffYOY_Shelikof <-> Sum_NearYOY_Kodiak,0,OffYOY_to_NearYOY
  #Spr_SST -> Sum_OffYOY_Shelikof,0, SST_to_OffYOY
  Sum_OffYOY_Shelikof -> recdevs, 1, offYOY_to_R
  Spr_SST -> recdevs, 1,SST_to_R
  Spr_SST -> Sum_OffYOY_Cond_Shelikof,0,SST_to_CondYOY
  Sum_OffYOY_Cond_Shelikof -> recdevs,1,CondYOY_to_R
  "

sem_iid = "
  # link, lag, param_name, start_value

  #internal data relation
  Spr_SST -> Spr_SST,1,AR_SST
  Wind_NS -> Wind_NS,1,AR_Wind
  Sum_LCopepod_Shelikof -> Sum_LCopepod_Shelikof,1, AR_Cope
  Sum_Euph_Kodiak -> Sum_Euph_Kodiak,1,AR_Euph
  Spr_Larvae_Shelikof - > Spr_Larvae_Shelikof,1,AR_larvae
  Sum_OffYOY_Cond_Shelikof -> Sum_OffYOY_Cond_Shelikof,1,AR_CondOffYOY
  Sum_OffYOY_Shelikof -> Sum_OffYOY_Shelikof,1,AR_OffYOY
  Sum_NearYOY_Kodiak -> Sum_NearYOY_Kodiak,1,AR_NearYOY
  Sum_Juv_EuphDiet -> Sum_Juv_EuphDiet,1,AR_JuvEuphDiet
  Fal_Adult_Cond_Fishery -> Fal_Adult_Cond_Fishery,1,AR_CondAd
  recdevs <-> recdevs, 0, sigmaR, 1
"

# first non fit
#devtools::install_github('James-Thorson-NOAA/dsem')
library(dsem)
family <- rep('normal',ncol(ESPdata_MSr))
DSEMcontrol <- dsem_control(use_REML=F, run_model=F,quiet=TRUE, getJointPrecision = TRUE,newton_loops=0)
fit_dsem = dsem(sem=sem_iid, tsdata=ts(ESPdata_MSr), family=family, control=DSEMcontrol )

# create new map & p
new_map = fit_dsem$tmb_inputs$map
new_map$lnsigma_j <- factor(rep(NA, length=length(new_map$lnsigma_j)))
new_pars <- fit_dsem$tmb_inputs$parameters
new_pars$lnsigma_j <- rep(log(0.1), length=length(new_pars$lnsigma_j))


library(TMB)
detach("package:Rceattle", unload = TRUE)
#devtools::build()
TMB::compile(file = 'src/TMB/ceattle_v01_11_dsem.cpp',
             PKG_CXXFLAGS = commandArgs(trailingOnly = TRUE),
             framework = "TMBad",
             safebounds = FALSE, safeunload = FALSE)
dyn.load(dynlib('src/TMB/ceattle_v01_11_dsem'))
library(Rceattle)

# The stock assessment inputs
map <- m0$obj$env$map
data <- m0$obj$env$data
pars <- m0$obj$env$parList(m0$opt$par)
random <- 'rec_dev'
# obj2 <- TMB::MakeADFun(data=data, parameters=pars, map=map,
#                        random=random, DLL='ceattle_v01_11',
#                        silent=TRUE)
# opt2 <- TMBhelper::fit_tmb(obj2)
# m0$opt$objective
# opt2$objective

# now build the joint DSEM model (SCEAM)
map2 <- c(map, new_map)
map2$rec_dev <- NULL
map2$R_ln_sd <- NULL
data2 <- c(data, fit_dsem$tmb_inputs$data)
pars2 <- c(pars,  fit_dsem$tmb_inputs$parameters)
pars2$rec_dev <- NULL
pars2$R_ln_sd <- NULL
random2 <- fit_dsem$tmb_inputs$random
stopifnot(ncol(pars$rec_dev)==nrow(pars2$x_tj))

obj_dsem <- TMB::MakeADFun(data=data2, parameters=pars2, map=map2,
                       random=random2, DLL='ceattle_v01_11_dsem',
                       silent=TRUE)
opt_dsem <- TMBhelper::fit_tmb(obj_dsem)
TMBhelper::check_estimability(obj_dsem)



