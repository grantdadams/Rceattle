data("GOApollock")
data_list <- GOApollock
data_list$Diet_comp_weights <- rep(1, data_list$nspp)

# data_list$R_sexr <- NULL
# data_list$est_sex_ratio <- NULL
# data_list$sex_ratio_sigma <- NULL
# data_list$aLW <- NULL
# data_list$maturity <- data_list$pmature
# data_list$pmature <- NULL
# data_list$weight <- data_list$wt
# data_list$wt <- NULL
#
# data_list$diet_data <- data_list$stom_prop_data
# data_list$stom_prop_data <- NULL
#
# ss_run <- Rceattle::fit_mod(data_list = data_list,
#                             inits = NULL, # Initial parameters = 0
#                             file = NULL, # Don't save
#                             estimateMode = 0, # Estimate
#                             random_rec = FALSE, # No random recruitment
#                             msmMode = 0, # Single species mode
#                             phase = TRUE,
#                             verbose = 1)
# plot_biomass(ss_run)

GOApollock <- data_list

usethis::use_data(GOApollock, overwrite = TRUE)
