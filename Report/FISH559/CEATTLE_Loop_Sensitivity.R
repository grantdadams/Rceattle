log_lik <- c()

load("~/Documents/GitHub/Rceattle/Report/FISH559/Models/ss_no_re.RData")
log_lik[1] <- mod_objects$opt$objective

for(i in 1:20){
  load(paste0("Report/FISH559/Models/ms_no_re_", i,".RData"))
  log_lik[i+1] <- mod_objects$opt$objective
}

loops <- c(0:20)

plot(x = loops, y = log_lik - min(log_lik), ylab = "Relative NLL", xlab = "M2 iterations", pch = 16)

# Sensitivity to other food
log_lik_of <- c()
other_food <- c()

for(i in 1:20){
  load(paste0("Report/FISH559/Models/Other Food/ms_no_re_5_optim_other_food_", i,".RData"))
  log_lik_of[i] <- mod_objects$opt$objective
  other_food[i] <- mod_objects$data_list$other_food[1]
}

loops <- c(1:20)
plot(x = other_food, y = log_lik_of - min(log_lik_of), ylab = "Relative NLL", xlab = "Other food", pch = 16)
