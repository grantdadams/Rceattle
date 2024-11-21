library(dplyr)
data("BS2017MS")
data_list <- BS2017MS
data_list$fleet_control <- data_list$fleet_control %>%
  rename(Estimate_index_sd = Estimate_survey_sd,
         Index_sd_prior = Survey_sd_prior)

BS2017MS <- data_list
usethis::use_data(BS2017MS, overwrite = T)
