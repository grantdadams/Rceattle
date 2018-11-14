#' Computes predation mortality using TMB or ADMB average N-at-age and ration estimates and TMB based diet data
#'
#' @param data_list List of data as used by TMB for CEATTLE
#' @param rep Report file including model outputs from TMB CEATTLE
#' @param tmp Report file including model outputs from ADMB CEATTLE
#' @param TMB Boolian indicating whether to use TMB (true) of ADMB (false) based average N-at-age and ration
#' @param admb_suit Boolian indicating whether or not to use suitability estimates from ADMB CEATTLE
#' @param admb_available_food Boolian indicating whether or not to use ADMB based available food estimates
#'
#' @return List of predation mortality based derived quantities
compare_pred_function <- function(data_list, rep = ms_rep, tmp = ms_tmp, TMB = TRUE, admb_suit = FALSE, admb_available_food = FALSE){

  # Step 1 - Extract ration and average N-at-age
  if(TMB == TRUE){
    AvgN <- list()
    AvgN[[1]] <- rep$AvgN[,,1]
    AvgN[[2]] <- rep$AvgN[,,2]
    AvgN[[3]] <- rep$AvgN[,,3]

    ration <- list()
    ration[[1]] <- rep$ration2Age[,,1]
    ration[[2]] <- rep$ration2Age[,,2]
    ration[[3]] <- rep$ration2Age[,,3]
  }
  if(TMB == FALSE){
    AvgN <- list()
    AvgN[[1]] <- tmp$AvgN_1
    AvgN[[2]] <- tmp$AvgN_2
    AvgN[[3]] <- tmp$AvgN_3

    ration <- list()
    ration[[1]] <- tmp$ration2_1
    ration[[2]] <- tmp$ration2_2
    ration[[3]] <- tmp$ration2_3
  }


  nyrs <- data_list$nyrs
  Uobs <- data_list$UobsAge


  # Step 2 - Get the percent other food in stomachs: (1 - sum(U)) / biomass of other food
  of_stomkir <- matrix(NA, 21, 3)
  for(pred in 1:3){
    for(pred_age in 1:data_list$nages[pred]){
      of_stomkir[pred_age, pred] <- (1 - sum( rep$stomKir[pred,,pred_age,,1], na.rm = T)) / data_list$other_food[1]
    }
  }


  # Stomach proportion: U/Biomass and sum(U/Biomass)
  stom_div_bio2 <- array(NA, dim = c(3,3,21,21, nyrs))
  suma_suit <- array(0, dim = c(nyrs, 21, 3))
  for(yr in 1:nyrs){
    for(prey in 1:3){
      for(pred in 1:3){
        for(prey_age in 1:data_list$nages[prey]){
          for(pred_age in 1:data_list$nages[pred]){
            stom_div_bio2[pred,prey,pred_age,prey_age, yr] <- rep$stomKir[pred,prey,pred_age,prey_age,yr] / ( AvgN[[prey]][yr, prey_age] * data_list$wt[yr, prey_age, prey] )
            suma_suit[yr,pred_age,pred] <- suma_suit[yr,pred_age,pred] + stom_div_bio2[pred,prey,pred_age,prey_age,yr]
          }
        }
      }
    }
  }

  # Step 3 - Suitability: (U/Biomass) / (sum(U/Biomass) + (1 - sum(U)) / biomass of other food)
  suit_main <- array(0, dim = c(3,3,21,21))
  suit_other <- matrix(0, 21, 3)
  for(pred in 1:3){
    for(pred_age in 1:data_list$nages[pred]){

      suit_other[pred_age, pred] = 1

      for(prey in 1:3){
        for(prey_age in 1:data_list$nages[prey]){
          for(yr in 1:nyrs){
            suit_main[pred,prey,pred_age,prey_age] <- suit_main[pred,prey,pred_age,prey_age] + (stom_div_bio2[pred,prey,pred_age,prey_age, yr] / (suma_suit[yr,pred_age,pred] + of_stomkir[pred_age,pred]))
          }
          suit_main[pred,prey,pred_age,prey_age] = suit_main[pred,prey,pred_age,prey_age] / nyrs
          suit_other[pred_age, pred] <- suit_other[pred_age, pred] - suit_main[pred,prey,pred_age,prey_age]
        }
      }
    }
  }


  # ADMB suitability
  suit_main_admb <- array(0, dim = c(3,3,21,21))
  suit_other_admb <- matrix(0, 21, 3)
  for(pred in 1:3){
    suit_other_admb[1:length( tmp[[paste("suit_other",pred, sep = "_")]] ), pred] <- replace(suit_other_admb[1:length( tmp[[paste("suit_other",pred, sep = "_")]] ), pred], values = tmp[[paste("suit_other",pred, sep = "_")]])
    for(pred_age in 1:data_list$nages[pred]){
      suit_main_admb[pred,,pred_age,] <- tmp[[paste("suit_main",pred, pred_age, sep = "_")]]
    }
  }

  if(admb_suit){
    suit_main <- suit_main_admb
    suit_other <- suit_other_admb
  }


  # Step 4 - Available food: (sum(suitability * biomass) + biomass of other prey * (1 - sum(suitability))) / nyrs
  avail_food <- array(0, dim = c(nyrs, 21, 3))
  for(pred in 1:3){
    for(pred_age in 1:data_list$nages[pred]){
      for(yr in 1:nyrs){
        tmp_othersuit = 0
        for(prey in 1:3){
          for(prey_age in 1:data_list$nages[prey]){

            avail_food[yr, pred_age, pred]  = avail_food[yr, pred_age, pred] + suit_main[pred,prey,pred_age,prey_age] * AvgN[[prey]][yr,prey_age] * data_list$wt[yr, prey_age, prey]
            tmp_othersuit = tmp_othersuit + suit_main[pred,prey,pred_age,prey_age]
          }
        }
        avail_food[yr, pred_age, pred] = avail_food[yr, pred_age, pred]  + data_list$other_food[1]*(1 - (tmp_othersuit))
      }
    }
  }


  if(admb_available_food){
    for(pred in 1:3){
      for(pred_age in 1:data_list$nages[pred]){
        for(yr in 1:nyrs){
          avail_food[yr, pred_age, pred] = tmp[[paste("avail_food",pred, sep = "_")]][yr, pred_age]
        }
      }
    }
  }


  # Step 5 -  Calculate M2: N * suitability * ration / available food
  M2 <- array(0, dim = c(nyrs, 21, 3))
  B_eaten <- array(0, dim = c(nyrs, 21, 3))
  for(prey in 1:3){
    for(prey_age in 1:data_list$nages[prey]){
      for(yr in 1:nyrs){
        for(pred in 1:3){
          for(pred_age in 1:data_list$nages[pred]){
            M2[yr, prey_age, prey] = M2[yr, prey_age, prey] + ((AvgN[[pred]][yr,pred_age] * ration[[pred]][yr,pred_age] * suit_main[pred,prey,pred_age,prey_age]) / avail_food[yr, pred_age, pred])
            B_eaten[yr, prey_age, prey] = B_eaten[yr, prey_age, prey] + (AvgN[[pred]][yr,pred_age] * ration[[pred]][yr,pred_age] * suit_main[pred,prey,pred_age,prey_age])
          }
        }
      }
    }
  }

  # Return results
  results <- list(R_M2 = M2, TMB_M2 = rep$M2, avail_food = avail_food, suit_main = suit_main, suit_other = suit_other, of_stomkir = of_stomkir, stom_div_bio2 = stom_div_bio2)
  return(results)
}


load("predation_data_and_outpurs.RData")

# Test 0 - All TMB
results <- compare_pred_function(data_list_ms , ms_rep, ms_tmp, TMB = TRUE, admb_suit = FALSE, admb_available_food = FALSE)
# Species 1
results$R_M2[1,1:12,1] # R
ms_tmp$M2_1[1,] # ADMB
ms_rep$M2[1,1:12,1] # TMB

sum(ms_rep$M2 - results$R_M2)


# Test 1 - ADMB average N-at-age, suitability, available food
results <- compare_pred_function(data_list_ms , ms_rep, ms_tmp, TMB = FALSE, admb_suit = TRUE, admb_available_food = TRUE)

# Species 1
results$R_M2[1,1:12,1] # R
ms_tmp$M2_1[1,] # ADMB
ms_rep$M2[1,1:12,1] # TMB


# Test 2 - ADMB average N-at-age and suitability, TMB based available food
results <- compare_pred_function(data_list_ms , ms_rep, ms_tmp, TMB = FALSE, admb_suit = TRUE, admb_available_food = FALSE)

# Species 1
results$R_M2[1,1:12,1] # R
ms_tmp$M2_1[1,] # ADMB
ms_rep$M2[1,1:12,1] # TMB


# Test 3 - ADMB average N-at-age - TMB based suitability and available food
results <- compare_pred_function(data_list_ms , ms_rep, ms_tmp, TMB = FALSE, admb_suit = FALSE, admb_available_food = FALSE)

# Species 1
results$R_M2[1,1:12,1] # R
ms_tmp$M2_1[1,] # ADMB
ms_rep$M2[1,1:12,1] # TMB


