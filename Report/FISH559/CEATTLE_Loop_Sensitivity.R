log_lik <- c()
log_lik_re <- c()

# Penalized like
load("Report/FISH559/Models/Optim 3/ss_no_re.RData")
log_lik[1] <- mod_objects$opt$objective

for(i in 1:20){
  load(paste0("Report/FISH559/Models/Optim 3/ms_no_re_", i,".RData"))
  log_lik[i+1] <- mod_objects$opt$objective
}

loops <- c(0:20)

# Random effects
load("Report/FISH559/Models/Optim 3/ss_re.RData")
log_lik_re[1] <- mod_objects$opt$objective

for(i in 1:10){
  load(paste0("Report/FISH559/Models/Optim 3/ms_re_", i,".RData"))
  log_lik_re[i+1] <- mod_objects$opt$objective
}

loops_re <- c(0:10)

# 5 optim
log_lik5 <- c()
load("Report/FISH559/Models/Optim 3/ss_no_re.RData")
log_lik5[1] <- mod_objects$opt$objective

for(i in 1:18){
  load(paste0("Report/FISH559/Models/Optim 5/ms_no_re_5optim_", i,".RData"))
  log_lik5[i+1] <- mod_objects$opt$objective
}

loops5 <- c(0:18)


# ADMB
log_lik_admb <- c()
loops_admb <- c(0,1,3,5,10,20)

dat <- as.character(read.delim(paste0("data/BS_SS_Files/ceattle.par"), header = FALSE)[1,])
dat <- as.numeric(strsplit(dat, split = " ")[[1]][11])
log_lik_admb[1] <- dat

for(i in 2:length(loops_admb)){
  dat <- as.character(read.delim(paste0("data/BS_MS_", loops_admb[i], "_Loops_Files/ceattle.par"), header = FALSE)[1,])
  dat <- as.numeric(strsplit(dat, split = " ")[[1]][11])
  log_lik_admb[i] <- dat
}



# Plot

plot(x = loops, y = log_lik - min(c(log_lik_admb, log_lik, log_lik_re)), ylab = "Relative NLL", xlab = "M2 iterations", type = "l", lwd = 3, ylim = c(0, max(c(log_lik_admb, log_lik, log_lik_re) - min(c(log_lik_admb, log_lik, log_lik_re)))))
lines(x = loops5, y = log_lik5 - min(c(log_lik, log_lik_admb, log_lik_re)), col = 1, pch = 16, lwd = 3)
lines(x = loops_re + 0.1, y = log_lik_re - min(c(log_lik, log_lik_admb, log_lik_re)), col = "grey60", pch = 16, lwd = 3)
lines(x = loops_admb + 0.2, y = log_lik_admb - min(c(log_lik, log_lik_admb, log_lik_re)), col = 3, pch = 16, lwd = 3, lty = 2)
legend("topright", legend = c("ADMB", "TMB", "TMB re"), col = c(3, 1, "grey60"), pch = rep(16, 3), bty = "n")











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
