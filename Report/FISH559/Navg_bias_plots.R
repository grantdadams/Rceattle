# Load models
load("Report/FISH559/Models/Optim 3/ms_no_re_10.RData")
ms_no_re <- mod_objects
ms_no_re$opt$AIC

load("Report/FISH559/Models/Optim 3/ms_re_10.RData")
ms_re <- mod_objects
ms_re$opt$AIC

load("Report/FISH559/Models/ms_no_re_n10.RData")
ms_n <- mod_objects
ms_n$opt$AIC



# Plot multispecies
load("C:/Users/Grant Adams/Documents/GitHub/Rceattle/data/BS_MS_10_Loops_Files_Corrected/CEATTLE_results.Rdata")
plot_biomass(ceattle_list = list(ms_no_re, ms_n),
             , file_name = "Report/FISH559/ms_avgn_bias",
             model_names = c("TMB pl Avg N", "TMB pl N"),
             line_col = c("#9B1D20", "#009FB7", adjustcolor( "red", alpha.f = 0.4)),
             species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder", rep(NA, 1)), lwd = 4)


plot_recruitment(ceattle_list = list(ms_no_re, ms_n),
                 file_name = "Report/FISH559/ms_avgn_bias",
                 model_names = c("TMB pl Avg N", "TMB pl N"),
                 line_col = c("#9B1D20", "#009FB7", adjustcolor( "red", alpha.f = .4)),
                 species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
                 ci_col = c("#F2BF79", "#ADCCBA", rep(adjustcolor( "red", alpha.f = 0), 1)), lwd = 4)



ms_no_re$quantities$M2[,1,1]
ms_n$quantities$M2[,1,1]


ms_no_re$quantities$AvgN[,4,1]
ms_n$quantities$AvgN[,4,1]

signif(ms_no_re$quantities$suit_main[1,1,1:12,1:12],6)
signif(ms_n$quantities$suit_main[1,1,1:12,1:12],6)


signif(ms_no_re$quantities$avail_food[1,,],6)
signif(ms_n$quantities$avail_food[1,,],6)
