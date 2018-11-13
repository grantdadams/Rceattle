# Compares predation functions using TMB or ADMB average age estimates and using TMB based diet and ration calculations
compare_pred_function <- function(data_list, rep = ms_rep, tmp = ms_tmp, TMB = TRUE, admb_suit = F){

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

  of_stomkir <- matrix(NA, 21, 3)
  for(pred in 1:3){
    for(pred_age in 1:21){
      of_stomkir[pred_age, pred] <- 1 - sum( rep$stomKir[pred,,pred_age,,1], na.rm = T)
    }
  }
  of_stomkir <- of_stomkir/data_list$other_food[1]

  # of_stomKir food should be good

  # Stomach proportion
  stom_div_bio2 <- array(NA, dim = c(3,3,21,21, nyrs))
  suma_suit <- array(0, dim = c(nyrs, 21, 3))
  for(yr in 1:nyrs){
    for(prey in 1:3){
      for(pred in 1:3){
        for(prey_age in 1:data_list$nages[prey]){
          for(pred_age in 1:data_list$nages[pred]){
            stom_div_bio2[pred,prey,pred_age,prey_age, yr] <- rep$stomKir[pred,prey,pred_age,prey_age,1] / (AvgN[[prey]][yr, prey_age] * data_list$wt[yr, prey_age, prey])
            suma_suit[yr,pred_age,pred] <- suma_suit[yr,pred_age,pred] + stom_div_bio2[pred,prey,pred_age,prey_age,yr]
          }
        }
      }
    }
  }
  suma_suit[1,,1]

  # Suitability
  suit_main <- array(0, dim = c(3,3,21,21))
  suit_other <- matrix(0, 21, 3)
  for(pred in 1:3){
    for(pred_age in 1:data_list$nages[pred]){
      suit_other[pred_age, pred] = 1
      for(prey in 1:3){
        for(prey_age in 1:data_list$nages[prey]){

          for(yr in 1:nyrs){
            suit_main[pred,prey,pred_age,prey_age] <- suit_main[pred,prey,pred_age,prey_age] + stom_div_bio2[pred,prey,pred_age,prey_age, yr] / (suma_suit[yr, pred_age, pred] + of_stomkir[pred_age,pred])
          }
          suit_main[pred,prey,pred_age,prey_age] = suit_main[pred,prey,pred_age,prey_age] / nyrs
          suit_other[pred_age, pred] <- suit_other[pred_age, pred] - suit_main[pred,prey,pred_age,prey_age]
        }
      }
    }
  }


  if(admb_suit){


    for(pred in 1:3){
      suit_other[1:length( tmp[[paste("suit_other",pred, sep = "_")]] ), pred] <- replace(suit_other[1:length( tmp[[paste("suit_other",pred, sep = "_")]] ), pred], values = tmp[[paste("suit_other",pred, sep = "_")]])

      for(pred_age in 1:data_list$nages[pred]){
        suit_main[pred,,pred_age,] <- tmp[[paste("suit_main",pred, pred_age, sep = "_")]]
      }
    }
  }


  # Compare
  tmp$suit_other_1
  tmp$suit_other_2
  tmp$suit_other_3

  # Available food
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

  tmp$avail_food_1[1,]
  tmp$avail_food_2[1,]
  tmp$avail_food_3[1,]

  avail_food[1,,1]
  avail_food[1,,2]
  avail_food[1,,3]

  # Calculate M2
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

  results <- list(R_M2 = M2, TMB_M2 = rep$M2)
  return(results)
}

results <- compare_pred_function(ms_run$data_list , ms_rep, ms_tmp, TMB = FALSE)
rounder = 4
# Species 1
results$R_M2[1,1:12,1]
ms_tmp$M2_1[1,]
ms_rep$M2[1,1:12,1]
