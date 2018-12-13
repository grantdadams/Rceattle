# Load models
load("Report/GSS_SS_admb.RData")
ss_admb <- mod_objects
ss_admb$opt$AIC

load("Report/GSS_SS_no_re.RData")
ss_no_re <- mod_objects
ss_no_re$opt$AIC

load("Report/GSS_SS_re.RData")
ss_re <- mod_objects
ss_re$opt$AIC

load("Report/GSS_MS_admb.RData")
ms_admb <- mod_objects
sum(ms_admb$quantities$jnll_comp) * 2 + 2 * (length(unlist(ms_admb$initial_params)) - 1)

load("Report/GSS_MS_no_re.RData")
ms_no_re <- mod_objects
ms_no_re$opt$AIC

load("Report/GSS_MS_re.RData")
ms_re <- mod_objects
ms_re$opt$AIC

aic_vec <- c(sum(ss_admb$quantities$jnll_comp) * 2 + 2 * (length(unlist(ss_admb$initial_params)) - 1), ss_no_re$opt$AIC, ss_re$opt$AIC, sum(ms_admb$quantities$jnll_comp) * 2 + 2 * (length(unlist(ms_admb$initial_params)) - 1), ms_no_re$opt$AIC, ms_re$opt$AIC)

sigma_r <- rbind( rep(0.707, 3) ,rep(0.707, 3) , ss_re$sdrep$value[which(names(ss_re$sdrep$value) == "r_sigma")], rep(0.707, 3), rep(0.707, 3), ms_re$sdrep$value[which(names(ms_re$sdrep$value) == "r_sigma")] )

sigma_r_sd <- rbind( rep(0, 3) ,rep(0, 3) , ss_re$sdrep$sd[which(names(ss_re$sdrep$value) == "r_sigma")], rep(0, 3), rep(0, 3), ms_re$sdrep$sd[which(names(ms_re$sdrep$value) == "r_sigma")] )

results <- cbind(aic_vec, sigma_r, sigma_r_sd)


results <- data.frame( AIC = aic_vec, sigmaR = )
ms_re$opt$AIC


# Plot biomass
plot_biomass(ceattle_list = list(ss_admb, ss_no_re, ss_re, ms_admb, ms_no_re, ms_re)
             , file_name = "Report/GSS/ss_runs", model_names = c("ADMB single spp", "TMB single spp", "TMB single spp re", rep(NA, 3)),
             line_col = c("#272727", "#9B1D20", "#009FB7", rep(NA, 3)),
             species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", rep(NA, 3)), lwd = 4)

plot_recruitment(ceattle_list = list(ss_admb, ss_no_re, ss_re, ms_admb, ms_no_re, ms_re)
                 , file_name = "Report/GSS/ss_runs", model_names = c("ADMB single spp", "TMB single spp", "TMB single spp re", rep(NA, 3)),
                 line_col = c("#272727", "#9B1D20", "#009FB7", rep(NA, 3)),
                 species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", rep(NA, 3)),
                 ci_col = c( "#7F7F7F", "#9B6F70", "#5BABB7", rep(adjustcolor( "red", alpha.f = 0), 3)), lwd = 4)


# Plot multispecies
plot_biomass(ceattle_list = list(ms_admb, ms_no_re, ms_re, ss_admb, ss_no_re, ss_re)
               , file_name = "Report/GSS/ms_runs",
             model_names = c("ADMB multi-spp", "TMB multi-spp", "TMB multi-spp re", rep(NA, 3)),
             line_col = c("#143642", "#F0C808", "#84BC9C", rep(NA, 3)),
             species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", rep(NA, 3)), lwd = 4)

plot_recruitment(ceattle_list = list(ms_admb, ms_no_re, ms_re, ss_admb, ss_no_re, ss_re)
                 , file_name = "Report/GSS/ms_runs",
                 model_names = c("ADMB multi-spp", "TMB multi-spp", "TMB multi-spp re", rep(NA, 3)),
                 line_col = c("#143642", "#F0C808", "#84BC9C", rep(NA, 3)),
                 species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", rep(NA, 3)),
                 ci_col = c( "#7F7F7F", "#F2BF79", "#ADCCBA", rep(adjustcolor( "red", alpha.f = 0), 3)), lwd = 4)


  # Plot RE
plot_biomass(ceattle_list = list(ms_re, ss_re, ms_admb, ms_no_re, ss_admb, ss_no_re)
             , file_name = "Report/GSS/re_runs",
             model_names = c("TMB multi-spp re", "TMB single-spp re", rep(NA, 4)),
             line_col = c("#84BC9C", "#009FB7", rep(NA, 4)),
             species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", rep(NA, 3)), lwd = 4)

plot_recruitment(ceattle_list = list(ms_re, ss_re, ms_admb, ms_no_re, ss_admb, ss_no_re)
                 , file_name = "Report/GSS/re_runs",
                 model_names = c("TMB multi-spp re", "TMB single-spp re", rep(NA, 4)),
                 line_col =  c("#84BC9C", "#009FB7", rep(NA, 4)),
                 species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", rep(NA, 3)),
                 ci_col = c( "#ADCCBA", "#5BABB7", rep(adjustcolor( "red", alpha.f = 0), 4)), lwd = 4)
