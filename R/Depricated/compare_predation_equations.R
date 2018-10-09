# Year 1
nyrs <- data_list$nyrs
Uobs <- data_list$UobsAge

of_stomkir <- matrix(NA, 21, 3)
for(pred in 1:3){
  for(pred_age in 1:21){
    of_stomkir[pred_age, pred] <- 1 - sum(Uobs[pred,,pred_age,],na.rm = T)
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
          stom_div_bio2[pred,prey,pred_age,prey_age, yr] <- Uobs[pred,prey,pred_age,prey_age] / (tmp[[paste0("AvgN_", prey)]][yr, prey_age] * data_list$wt[yr, prey_age, prey])
          suma_suit[yr,pred_age,pred] <- suma_suit[yr,pred_age,pred] + stom_div_bio2[pred,prey,pred_age,prey_age,yr]
        }
      }
    }
  }
}

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

suit_other

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

          avail_food[yr, pred_age, pred]  = avail_food[yr, pred_age, pred] + suit_main[pred,prey,pred_age,prey_age] * tmp[[paste0("AvgN_", prey)]][yr,prey_age] * data_list$wt[yr, prey_age, prey]
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

# Calculate M2
M2 <- array(0, dim = c(nyrs, 21, 3))
B_eaten <- array(0, dim = c(nyrs, 21, 3))
for(prey in 1:3){
  for(prey_age in 1:data_list$nages[prey]){
    for(yr in 1:nyrs){
      for(pred in 1:3){
        for(pred_age in 1:data_list$nages[pred]){
          M2[yr, prey_age, prey] = M2[yr, prey_age, prey] + ((tmp[[paste0("AvgN_", pred)]][yr,pred_age] * tmp[[paste0("ration2_", pred)]][yr,pred_age]*suit_main[pred,prey,pred_age,prey_age])/avail_food[yr, pred_age, pred])
          B_eaten[yr, prey_age, prey] = B_eaten[yr, prey_age, prey] + (tmp[[paste0("AvgN_", pred)]][yr,pred_age] * tmp[[paste0("ration2_", pred)]][yr,pred_age]*suit_main[pred,prey,pred_age,prey_age])
        }
      }
    }
  }
}

