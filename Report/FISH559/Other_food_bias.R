
library(Rceattle)
mod_list <- list()
other_food_list <- c()
for(i in 1:40){
  load(paste0("Report/FISH559/Models/Other Food/ms_no_re_5_optim_other_food_",i,".RData"))
  print(i)
  mod_list[[i]] <- mod_objects
  other_food_list[i] <- mod_objects$data_list$other_food[1]
}

colfunc <- colorRampPalette(c("black", "red"))


other_food_percent <- seq(0.05, 2, by = 0.05 )
cols <- rev(colfunc(length(other_food_percent)))

plot_biomass(ceattle_list = mod_list,
             file_name = "Report/FISH559/other_food_b",
             model_names = paste0(other_food_percent * 100, "%"),
             line_col = cols,
             species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"), lwd = 3)


plot_recruitment(ceattle_list = mod_list,
                 file_name = "Report/FISH559/other_food_b",
                 model_names = paste0(other_food_percent * 100, "%"),
                 line_col = cols,
                 species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
                 ci_col = NULL, lwd = 3)


like <- sapply(mod_list, function(x) x$opt$objective)

plot( x = paste0(other_food_percent * 100), y = like - min(like), xlab = "% other food", ylab = "Relative NLL", pch = 16)

