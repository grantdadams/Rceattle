load("data/BS_MS_3_Loops_Files/CEATTLE_results.Rdata")
pollock_suitability_avgN <- tmp$AvgNsuit_1
pollock_M2_avgN <- tmp$AvgNM2_1

colnames(pollock_M2_avgN) <- paste0("Age_", 1:ncol(pollock_M2_avgN))
colnames(pollock_suitability_avgN) <- paste0("Age_", 1:ncol(pollock_suitability_avgN))

pollock_suitability_avgN[1:10,]
pollock_M2_avgN[1:10,]
