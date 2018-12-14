load("data/BS_SS_Files/2017_assessment_data_list.RData")
data_list_ss$msmMode=1

load("data/BS_MS_Files/2017_assessment_data_list.RData")


length(data_list_ms)
length(data_list_ss)
same <- list()


for(i in 1:length(data_list_ms)){
  same[[i]] <- data_list_ms[[i]] == data_list_ss[[i]]
  names(same)[i] <- names(data_list_ms)[i]
}


res <- lapply(same, function(ch) grep(FALSE, ch))
